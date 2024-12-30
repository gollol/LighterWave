#include <fstream>
#include <functional>
#include <lightwave.hpp>
#include <stack>

namespace lightwave {

static void addToAtomicFloat(std::atomic<float> &var, float val) {
    auto current = var.load();
    while (!var.compare_exchange_weak(current, current + val))
        ;
}

inline float logistic(float x) { return 1 / (1 + std::exp(-x)); }

// Implements the stochastic-gradient-based Adam optimizer [Kingma and Ba 2014]
class AdamOptimizer {
public:
    AdamOptimizer(float learningRate, int batchSize = 1, float epsilon = 1e-08f,
                  float beta1 = 0.9f, float beta2 = 0.999f) {
        m_hparams = { learningRate, batchSize, epsilon, beta1, beta2 };
    }

    AdamOptimizer &operator=(const AdamOptimizer &arg) {
        m_state   = arg.m_state;
        m_hparams = arg.m_hparams;
        return *this;
    }

    AdamOptimizer(const AdamOptimizer &arg) { *this = arg; }

    void append(float gradient, float statisticalWeight) {
        m_state.batchGradient += gradient * statisticalWeight;
        m_state.batchAccumulation += statisticalWeight;

        if (m_state.batchAccumulation > m_hparams.batchSize) {

            step(m_state.batchGradient / m_state.batchAccumulation);

            m_state.batchGradient     = 0;
            m_state.batchAccumulation = 0;
        }
    }

    void step(float gradient) {
        ++m_state.iter;
        m_state.prevsecondmoment = m_state.secondMoment;

        float actualLearningRate =
            m_hparams.learningRate *
            safe_sqrt(1 - std::pow(m_hparams.beta2, m_state.iter)) /
            (1 - std::pow(m_hparams.beta1, m_state.iter));
        m_state.firstMoment = m_hparams.beta1 * m_state.firstMoment +
                              (1 - m_hparams.beta1) * gradient;
        m_state.secondMoment = m_hparams.beta2 * m_state.secondMoment +
                               (1 - m_hparams.beta2) * gradient * gradient;
        m_state.variable -=
            actualLearningRate * m_state.firstMoment /
            (safe_sqrt(m_state.secondMoment) + m_hparams.epsilon);

        // Clamp the variable to the range [-20, 20] as a safeguard to avoid
        // numerical instability: since the sigmoid involves the exponential of
        // the variable, value of -20 or 20 already yield in *extremely* small
        // and large results that are pretty much never necessary in practice.
        m_state.variable = std::min(std::max(m_state.variable, -20.0f), 20.0f);
    }

    float variable() const { return m_state.variable; }

private:
    struct State {
        int iter               = 0;
        float firstMoment      = 0;
        float secondMoment     = 0;
        float prevsecondmoment = secondMoment;
        float variable         = 0;

        float batchAccumulation = 0;
        float batchGradient     = 0;
    } m_state;

    struct Hyperparameters {
        float learningRate;
        int batchSize;
        float epsilon;
        float beta1;
        float beta2;
    } m_hparams;
};

enum ENee {
    ENever,
    EKickstart,
    EAlways,
};

enum class EBsdfSamplingFractionLoss {
    ENone,
    EKL,
    EVariance,
};

enum class ESpatialFilter {
    ENearest,
    EStochasticBox,
    EBox,
};

enum class EDirectionalFilter {
    ENearest,
    EBox,
};

class QuadTreeNode {
public:
    QuadTreeNode() {
        m_children = {};
        for (size_t i = 0; i < m_sum.size(); ++i) {
            m_sum[i].store(0, std::memory_order_relaxed);
        }
    }

    void setSum(int index, float val) {
        m_sum[index].store(val, std::memory_order_relaxed);
    }

    float sum(int index) const {
        return m_sum[index].load(std::memory_order_relaxed);
    }

    void copyFrom(const QuadTreeNode &arg) {
        for (int i = 0; i < 4; ++i) {
            setSum(i, arg.sum(i));
            m_children[i] = arg.m_children[i];
        }
    }

    QuadTreeNode(const QuadTreeNode &arg) { copyFrom(arg); }

    QuadTreeNode &operator=(const QuadTreeNode &arg) {
        copyFrom(arg);
        return *this;
    }

    void setChild(int idx, uint16_t val) { m_children[idx] = val; }

    uint16_t child(int idx) const { return m_children[idx]; }

    void setSum(float val) {
        for (int i = 0; i < 4; ++i) {
            setSum(i, val);
        }
    }

    int childIndex(Point2 &p) const {
        int res = 0;
        for (int i = 0; i < 2; ++i) {
            if (p[i] < 0.5f) {
                p[i] *= 2;
            } else {
                p[i] = (p[i] - 0.5f) * 2;
                res |= 1 << i;
            }
        }

        return res;
    }

    // Evaluates the directional irradiance *sum density* (i.e. sum / area) at a
    // given location p. To obtain radiance, the sum density (result of this
    // function) must be divided by the total statistical weight of the
    // estimates that were summed up.
    float eval(Point2 &p, const std::vector<QuadTreeNode> &nodes) const {
        assert(p.x() >= 0 && p.x() <= 1 && p.y() >= 0 && p.y() <= 1);
        const int index = childIndex(p);
        if (isLeaf(index)) {
            return 4 * sum(index);
        } else {
            return 4 * nodes[child(index)].eval(p, nodes);
        }
    }

    float pdf(Point2 &p, const std::vector<QuadTreeNode> &nodes) const {
        assert(p.x() >= 0 && p.x() <= 1 && p.y() >= 0 && p.y() <= 1);
        const int index = childIndex(p);
        if (!(sum(index) > 0)) {
            return 0;
        }

        const float factor =
            4 * sum(index) / (sum(0) + sum(1) + sum(2) + sum(3));
        if (isLeaf(index)) {
            return factor;
        } else {
            return factor * nodes[child(index)].pdf(p, nodes);
        }
    }

    int depthAt(Point2 &p, const std::vector<QuadTreeNode> &nodes) const {
        assert(p.x() >= 0 && p.x() <= 1 && p.y() >= 0 && p.y() <= 1);
        const int index = childIndex(p);
        if (isLeaf(index)) {
            return 1;
        } else {
            return 1 + nodes[child(index)].depthAt(p, nodes);
        }
    }

    Point2 sample(Sampler *sampler,
                  const std::vector<QuadTreeNode> &nodes) const {
        int index = 0;

        float topLeft  = sum(0);
        float topRight = sum(1);
        float partial  = topLeft + sum(2);
        float total    = partial + topRight + sum(3);

        // Should only happen when there are numerical instabilities.
        if (!(total > 0.0f)) {
            return sampler->next2D();
        }

        float boundary = partial / total;
        Point2 origin  = Point2{ 0.0f, 0.0f };

        float sample = sampler->next();

        if (sample < boundary) {
            assert(partial > 0);
            sample /= boundary;
            boundary = topLeft / partial;
        } else {
            partial = total - partial;
            assert(partial > 0);
            origin.x() = 0.5f;
            sample     = (sample - boundary) / (1.0f - boundary);
            boundary   = topRight / partial;
            index |= 1 << 0;
        }

        if (sample < boundary) {
            sample /= boundary;
        } else {
            origin.y() = 0.5f;
            sample     = (sample - boundary) / (1.0f - boundary);
            index |= 1 << 1;
        }

        if (isLeaf(index)) {
            Point2 sample = sampler->next2D();
            return Point2(origin.x() + 0.5f * sample.x(),
                          origin.y() + 0.5f * sample.y());
        } else {
            Point2 sample = nodes[child(index)].sample(sampler, nodes);
            return Point2(origin.x() + 0.5f * sample.x(),
                          origin.y() + 0.5f * sample.y());
        }
    }

    void record(Point2 &p, float irradiance, std::vector<QuadTreeNode> &nodes) {
        assert(p.x() >= 0 && p.x() <= 1 && p.y() >= 0 && p.y() <= 1);
        int index = childIndex(p);

        if (isLeaf(index)) {
            addToAtomicFloat(m_sum[index], irradiance);
        } else {
            nodes[child(index)].record(p, irradiance, nodes);
        }
    }

    float computeOverlappingArea(const Point2 &min1, const Point2 &max1,
                                 const Point2 &min2, const Point2 &max2) {
        float lengths[2];

        float size = max2[0] - min2[0];
        if (max2[0] + size < max1[0]) {
            lengths[0] = std::max(std::min(max1[0], max2[0] - 2 * size) -
                                      std::max(min1[0], min2[0] - 2 * size),
                                  0.0f);

        } else if (min2[0] - size > min1[0]) {
            lengths[0] = std::max(std::min(max1[0], max2[0] + 2 * size) -
                                      std::max(min1[0], min2[0] + 2 * size),
                                  0.0f);

        } else {
            lengths[0] = std::max(
                std::min(max1[0], max2[0]) - std::max(min1[0], min2[0]), 0.0f);
        }

        lengths[1] = std::max(
            std::min(max1[1], max2[1]) - std::max(min1[1], min2[1]), 0.0f);
        return lengths[0] * lengths[1];
    }

    void record(const Point2 &origin, float size, Point2 nodeOrigin,
                float nodeSize, float value, std::vector<QuadTreeNode> &nodes) {
        float childSize = nodeSize / 2;
        for (int i = 0; i < 4; ++i) {
            Point2 childOrigin = nodeOrigin;
            if (i & 1) {
                childOrigin[0] += childSize;
            }
            if (i & 2) {
                childOrigin[1] += childSize;
            }

            float w = computeOverlappingArea(
                origin,
                Point2(origin.x() + size, origin.y() + size),
                childOrigin,
                Point2(childOrigin.x() + childSize,
                       childOrigin.y() + childSize));
            if (w > 0.0f) {
                if (isLeaf(i)) {
                    addToAtomicFloat(m_sum[i], value * w);
                } else {
                    nodes[child(i)].record(
                        origin, size, childOrigin, childSize, value, nodes);
                }
            }
        }
    }

    bool isLeaf(int index) const { return child(index) == 0; }

    // Ensure that each quadtree node's sum of irradiance estimates
    // equals that of all its children.
    void build(std::vector<QuadTreeNode> &nodes) {
        for (int i = 0; i < 4; ++i) {
            // During sampling, all irradiance estimates are accumulated in
            // the leaves, so the leaves are built by definition.
            if (isLeaf(i)) {
                continue;
            }

            QuadTreeNode &c = nodes[child(i)];

            // Recursively build each child such that their sum becomes valid...
            c.build(nodes);

            // ...then sum up the children's sums.
            float sum = 0;
            for (int j = 0; j < 4; ++j) {
                sum += c.sum(j);
            }
            setSum(i, sum);
        }
    }

private:
    std::array<std::atomic<float>, 4> m_sum;
    std::array<uint16_t, 4> m_children;
};

class DTree {
public:
    DTree() {
        m_atomic.sum.store(0, std::memory_order_relaxed);
        m_maxDepth = 0;
        m_nodes.emplace_back();
        m_nodes.front().setSum(0.0f);
    }

    const QuadTreeNode &node(size_t i) const { return m_nodes[i]; }

    float mean() const {
        if (m_atomic.statisticalWeight == 0) {
            return 0;
        }
        const float factor = 1 / (Pi * 4 * m_atomic.statisticalWeight);
        return factor * m_atomic.sum;
    }

    void recordIrradiance(Point2 p, float irradiance, float statisticalWeight,
                          EDirectionalFilter directionalFilter) {
        if (std::isfinite(statisticalWeight) && statisticalWeight > 0) {
            addToAtomicFloat(m_atomic.statisticalWeight, statisticalWeight);

            if (std::isfinite(irradiance) && irradiance > 0) {
                if (directionalFilter == EDirectionalFilter::ENearest) {
                    m_nodes[0].record(
                        p, irradiance * statisticalWeight, m_nodes);
                } else {
                    int depth  = depthAt(p);
                    float size = std::pow(0.5f, depth);

                    Point2 origin = p;
                    origin.x() -= size / 2;
                    origin.y() -= size / 2;
                    m_nodes[0].record(origin,
                                      size,
                                      Point2(0.0f),
                                      1.0f,
                                      irradiance * statisticalWeight /
                                          (size * size),
                                      m_nodes);
                }
            }
        }
    }

    float pdf(Point2 p) const {
        if (!(mean() > 0)) {
            return 1 / (4 * Pi);
        }

        return m_nodes[0].pdf(p, m_nodes) / (4 * Pi);
    }

    int depthAt(Point2 p) const { return m_nodes[0].depthAt(p, m_nodes); }

    int depth() const { return m_maxDepth; }

    Point2 sample(Sampler *sampler) const {
        if (!(mean() > 0)) {
            return sampler->next2D();
        }

        Point2 res = m_nodes[0].sample(sampler, m_nodes);

        res.x() = clamp(res.x(), 0.0f, 1.0f);
        res.y() = clamp(res.y(), 0.0f, 1.0f);

        return res;
    }

    size_t numNodes() const { return m_nodes.size(); }

    float statisticalWeight() const { return m_atomic.statisticalWeight; }

    void setStatisticalWeight(float statisticalWeight) {
        m_atomic.statisticalWeight = statisticalWeight;
    }

    void reset(const DTree &previousDTree, int newMaxDepth,
               float subdivisionThreshold) {
        m_atomic   = Atomic{};
        m_maxDepth = 0;
        m_nodes.clear();
        m_nodes.emplace_back();

        struct StackNode {
            size_t nodeIndex;
            size_t otherNodeIndex;
            const DTree *otherDTree;
            int depth;
        };

        std::stack<StackNode> nodeIndices;
        nodeIndices.push({ 0, 0, &previousDTree, 1 });

        const float total = previousDTree.m_atomic.sum;

        // Create the topology of the new DTree to be the refined version
        // of the previous DTree. Subdivision is recursive if enough energy is
        // there.
        while (!nodeIndices.empty()) {
            StackNode sNode = nodeIndices.top();
            nodeIndices.pop();

            m_maxDepth = std::max(m_maxDepth, sNode.depth);

            for (int i = 0; i < 4; ++i) {
                const QuadTreeNode &otherNode =
                    sNode.otherDTree->m_nodes[sNode.otherNodeIndex];
                const float fraction = total > 0 ? (otherNode.sum(i) / total)
                                                 : std::pow(0.25f, sNode.depth);
                assert(fraction <= 1.0f + Epsilon);

                if (sNode.depth < newMaxDepth &&
                    fraction > subdivisionThreshold) {
                    if (!otherNode.isLeaf(i)) {
                        assert(sNode.otherDTree == &previousDTree);
                        nodeIndices.push({ m_nodes.size(),
                                           otherNode.child(i),
                                           &previousDTree,
                                           sNode.depth + 1 });
                    } else {
                        nodeIndices.push({ m_nodes.size(),
                                           m_nodes.size(),
                                           this,
                                           sNode.depth + 1 });
                    }

                    m_nodes[sNode.nodeIndex].setChild(
                        i, static_cast<uint16_t>(m_nodes.size()));
                    m_nodes.emplace_back();
                    m_nodes.back().setSum(otherNode.sum(i) / 4);

                    if (m_nodes.size() > std::numeric_limits<uint16_t>::max()) {
                        logger(EWarn,
                               "DTreeWrapper hit maximum children count.");
                        nodeIndices = std::stack<StackNode>();
                        break;
                    }
                }
            }
        }

        // Uncomment once memory becomes an issue.
        // m_nodes.shrink_to_fit();

        for (auto &node : m_nodes) {
            node.setSum(0);
        }
    }

    size_t approxMemoryFootprint() const {
        return m_nodes.capacity() * sizeof(QuadTreeNode) + sizeof(*this);
    }

    void build() {
        auto &root = m_nodes[0];

        // Build the quadtree recursively, starting from its root.
        root.build(m_nodes);

        // Ensure that the overall sum of irradiance estimates equals
        // the sum of irradiance estimates found in the quadtree.
        float sum = 0;
        for (int i = 0; i < 4; ++i) {
            sum += root.sum(i);
        }
        m_atomic.sum.store(sum);
    }

private:
    std::vector<QuadTreeNode> m_nodes;

    struct Atomic {
        Atomic() {
            sum.store(0, std::memory_order_relaxed);
            statisticalWeight.store(0, std::memory_order_relaxed);
        }

        Atomic(const Atomic &arg) { *this = arg; }

        Atomic &operator=(const Atomic &arg) {
            sum.store(arg.sum.load(std::memory_order_relaxed),
                      std::memory_order_relaxed);
            statisticalWeight.store(
                arg.statisticalWeight.load(std::memory_order_relaxed),
                std::memory_order_relaxed);
            return *this;
        }

        std::atomic<float> sum;
        std::atomic<float> statisticalWeight;

    } m_atomic;

    int m_maxDepth;
};

struct DTreeRecord {
    Vector d;
    float radiance, product;
    float woPdf, bsdfPdf, dTreePdf;
    float statisticalWeight;
    bool isDelta;
};

struct DTreeWrapper {
public:
    DTreeWrapper() {}

    void record(const DTreeRecord &rec, EDirectionalFilter directionalFilter,
                EBsdfSamplingFractionLoss bsdfSamplingFractionLoss) {
        if (!rec.isDelta) {
            float irradiance = rec.radiance / rec.woPdf;
            building.recordIrradiance(dirToCanonical(rec.d),
                                      irradiance,
                                      rec.statisticalWeight,
                                      directionalFilter);
        }

        if (bsdfSamplingFractionLoss != EBsdfSamplingFractionLoss::ENone &&
            rec.product > 0) {
            optimizeBsdfSamplingFraction(rec,
                                         bsdfSamplingFractionLoss ==
                                                 EBsdfSamplingFractionLoss::EKL
                                             ? 1.0f
                                             : 2.0f);
        }
    }

    static Vector canonicalToDir(Point2 p) {
        const float cosTheta = 2 * p.x() - 1;
        const float phi      = 2 * Pi * p.y();

        const float sinTheta = safe_sqrt(1 - cosTheta * cosTheta);
        float sinPhi         = sinf(phi);
        float cosPhi         = cosf(phi);

        return { sinTheta * cosPhi, sinTheta * sinPhi, cosTheta };
    }

    static Point2 dirToCanonical(const Vector &d) {
        if (!std::isfinite(d.x()) || !std::isfinite(d.y()) ||
            !std::isfinite(d.z())) {
            return { 0, 0 };
        }

        const float cosTheta = std::min(std::max(d.z(), -1.0f), 1.0f);
        float phi            = std::atan2(d.y(), d.x());
        while (phi < 0)
            phi += 2.0 * Pi;

        return { (cosTheta + 1) / 2, phi / (2 * Pi) };
    }

    void build() {
        building.build();
        sampling = building;
    }

    void reset(int maxDepth, float subdivisionThreshold) {
        building.reset(sampling, maxDepth, subdivisionThreshold);
    }

    Vector sample(Sampler *sampler) const {
        return canonicalToDir(sampling.sample(sampler));
    }

    float pdf(const Vector &dir) const {
        return sampling.pdf(dirToCanonical(dir));
    }

    float diff(const DTreeWrapper &other) const { return 0.0f; }

    int depth() const { return sampling.depth(); }

    size_t numNodes() const { return sampling.numNodes(); }

    float meanRadiance() const { return sampling.mean(); }

    float statisticalWeight() const { return sampling.statisticalWeight(); }

    float statisticalWeightBuilding() const {
        return building.statisticalWeight();
    }

    void setStatisticalWeightBuilding(float statisticalWeight) {
        building.setStatisticalWeight(statisticalWeight);
    }

    size_t approxMemoryFootprint() const {
        return building.approxMemoryFootprint() +
               sampling.approxMemoryFootprint();
    }

    inline float bsdfSamplingFraction(float variable) const {
        return logistic(variable);
    }

    inline float dBsdfSamplingFraction_dVariable(float variable) const {
        float fraction = bsdfSamplingFraction(variable);
        return fraction * (1 - fraction);
    }

    inline float bsdfSamplingFraction() const {
        return bsdfSamplingFraction(bsdfSamplingFractionOptimizer.variable());
    }

    void optimizeBsdfSamplingFraction(const DTreeRecord &rec,
                                      float ratioPower) {
        m_lock.lock();

        // GRADIENT COMPUTATION
        float variable         = bsdfSamplingFractionOptimizer.variable();
        float samplingFraction = bsdfSamplingFraction(variable);

        // Loss gradient w.r.t. sampling fraction
        float mixPdf = samplingFraction * rec.bsdfPdf +
                       (1 - samplingFraction) * rec.dTreePdf;
        float ratio = std::pow(rec.product / mixPdf, ratioPower);
        float dLoss_dSamplingFraction =
            -ratio / rec.woPdf * (rec.bsdfPdf - rec.dTreePdf);

        // Chain rule to get loss gradient w.r.t. trainable variable
        float dLoss_dVariable =
            dLoss_dSamplingFraction * dBsdfSamplingFraction_dVariable(variable);

        // We want some regularization such that our parameter does not become
        // too big. We use l2 regularization, resulting in the following linear
        // gradient.
        float l2RegGradient = 0.01f * variable;

        float lossGradient = l2RegGradient + dLoss_dVariable;

        // ADAM GRADIENT DESCENT
        bsdfSamplingFractionOptimizer.append(lossGradient,
                                             rec.statisticalWeight);

        m_lock.unlock();
    }

    Color getLrEstimate() const {
        return weight > 0.0f ? firstMoment / weight : Color(0.0f);
    }

    void splatLrEstimate(const Color &sum, float wght) {
        atomicAdd(firstMoment, sum);
        atomicAdd(weight, wght);
    }

    void setLrEstimate(const Color sum, float wght) {
        firstMoment = sum;
        weight      = wght;
    }
    Color firstMoment = Color(0.0f);
    float weight      = 0.0f;

private:
    DTree building;
    DTree sampling;

    AdamOptimizer bsdfSamplingFractionOptimizer{ 0.01f };

    class SpinLock {
    public:
        SpinLock() { m_mutex.clear(std::memory_order_release); }

        SpinLock(const SpinLock &other) {
            m_mutex.clear(std::memory_order_release);
        }
        SpinLock &operator=(const SpinLock &other) { return *this; }

        void lock() {
            while (m_mutex.test_and_set(std::memory_order_acquire)) {
            }
        }

        void unlock() { m_mutex.clear(std::memory_order_release); }

    private:
        std::atomic_flag m_mutex;
    } m_lock;
};

struct STreeNode {
    STreeNode() {
        children = {};
        isLeaf   = true;
        axis     = 0;
    }

    int childIndex(Point &p) const {
        if (p[axis] < 0.5f) {
            p[axis] *= 2;
            return 0;
        } else {
            p[axis] = (p[axis] - 0.5f) * 2;
            return 1;
        }
    }

    int nodeIndex(Point &p) const { return children[childIndex(p)]; }

    DTreeWrapper *dTreeWrapper(Point &p, Vector &size,
                               std::vector<STreeNode> &nodes) {
        assert(p[axis] >= -Epsilon && p[axis] <= 1 + Epsilon);
        p[axis] = clamp(p[axis], 0.0f, 1.0f);
        if (isLeaf) {
            return &dTree;
        } else {
            size[axis] /= 2;
            return nodes[nodeIndex(p)].dTreeWrapper(p, size, nodes);
        }
    }

    const DTreeWrapper *dTreeWrapper() const { return &dTree; }

    int depth(Point &p, const std::vector<STreeNode> &nodes) const {
        assert(p[axis] >= 0 && p[axis] <= 1);
        if (isLeaf) {
            return 1;
        } else {
            return 1 + nodes[nodeIndex(p)].depth(p, nodes);
        }
    }

    int depth(const std::vector<STreeNode> &nodes) const {
        int result = 1;

        if (!isLeaf) {
            for (auto c : children) {
                result = std::max(result, 1 + nodes[c].depth(nodes));
            }
        }

        return result;
    }

    void forEachLeaf(
        std::function<void(const DTreeWrapper *, const Point &, const Vector &)>
            func,
        Point p, Vector size, const std::vector<STreeNode> &nodes) const {

        if (isLeaf) {
            func(&dTree, p, size);
        } else {
            size[axis] /= 2;
            for (int i = 0; i < 2; ++i) {
                Point childP = p;
                if (i == 1) {
                    childP[axis] += size[axis];
                }

                nodes[children[i]].forEachLeaf(func, childP, size, nodes);
            }
        }
    }

    float computeOverlappingVolume(const Point &min1, const Point &max1,
                                   const Point &min2, const Point &max2) {
        float lengths[3];
        for (int i = 0; i < 3; ++i) {
            lengths[i] = std::max(
                std::min(max1[i], max2[i]) - std::max(min1[i], min2[i]), 0.0f);
        }
        return lengths[0] * lengths[1] * lengths[2];
    }

    void record(const Point &min1, const Point &max1, Point min2, Vector size2,
                const DTreeRecord &rec, EDirectionalFilter directionalFilter,
                EBsdfSamplingFractionLoss bsdfSamplingFractionLoss,
                std::vector<STreeNode> &nodes) {
        float w = computeOverlappingVolume(min1, max1, min2, min2 + size2);
        if (w > 0) {
            if (isLeaf) {
                dTree.record({ rec.d,
                               rec.radiance,
                               rec.product,
                               rec.woPdf,
                               rec.bsdfPdf,
                               rec.dTreePdf,
                               rec.statisticalWeight * w,
                               rec.isDelta },
                             directionalFilter,
                             bsdfSamplingFractionLoss);
            } else {
                size2[axis] /= 2;
                for (int i = 0; i < 2; ++i) {
                    if (i & 1) {
                        min2[axis] += size2[axis];
                    }

                    nodes[children[i]].record(min1,
                                              max1,
                                              min2,
                                              size2,
                                              rec,
                                              directionalFilter,
                                              bsdfSamplingFractionLoss,
                                              nodes);
                }
            }
        }
    }

    bool isLeaf;
    DTreeWrapper dTree;
    int axis;
    std::array<uint32_t, 2> children;
};

class STree {
public:
    STree(const Bounds &aabb) {
        clear();
        logger(EInfo, "STree created.");

        // Enlarge AABB to turn it into a cube. This has the effect
        // of nicer hierarchical subdivisions.
        Vector size   = aabb.max() - aabb.min();
        float maxSize = std::max(std::max(size.x(), size.y()), size.z());
        m_aabb        = Bounds(aabb.min(), aabb.min() + Vector(maxSize));
    }

    void clear() {
        m_nodes.clear();
        m_nodes.emplace_back();
    }

    void subdivideAll() {
        int nNodes = (int) m_nodes.size();
        for (int i = 0; i < nNodes; ++i) {
            if (m_nodes[i].isLeaf) {
                subdivide(i, m_nodes);
            }
        }
    }

    void subdivide(int nodeIdx, std::vector<STreeNode> &nodes) {
        // Add 2 child nodes
        nodes.resize(nodes.size() + 2);

        if (nodes.size() > std::numeric_limits<uint32_t>::max()) {
            logger(EWarn, "DTreeWrapper hit maximum children count.");
            return;
        }

        STreeNode &cur = nodes[nodeIdx];
        for (int i = 0; i < 2; ++i) {
            uint32_t idx     = (uint32_t) nodes.size() - 2 + i;
            cur.children[i]  = idx;
            nodes[idx].axis  = (cur.axis + 1) % 3;
            nodes[idx].dTree = cur.dTree;
            nodes[idx].dTree.setStatisticalWeightBuilding(
                nodes[idx].dTree.statisticalWeightBuilding() / 2);
            // eyeball it //TODO: check
            nodes[idx].dTree.setLrEstimate(cur.dTree.firstMoment / 2,
                                           cur.dTree.weight / 2);
        }
        cur.isLeaf = false;
        cur.dTree  = {}; // Reset to an empty dtree to save memory.
    }

    DTreeWrapper *dTreeWrapper(Point p, Vector &size) {
        size = m_aabb.max() - m_aabb.min();
        p    = Point(p - m_aabb.min());
        p.x() /= size.x();
        p.y() /= size.y();
        p.z() /= size.z();
        return m_nodes[0].dTreeWrapper(p, size, m_nodes);
    }

    DTreeWrapper *dTreeWrapper(Point p) {
        Vector size;
        return dTreeWrapper(p, size);
    }

    void forEachDTreeWrapperConst(
        std::function<void(const DTreeWrapper *)> func) const {
        for (auto &node : m_nodes) {
            if (node.isLeaf) {
                func(&node.dTree);
            }
        }
    }

    void forEachDTreeWrapperConstP(
        std::function<void(const DTreeWrapper *, const Point &, const Vector &)>
            func) const {
        m_nodes[0].forEachLeaf(
            func, m_aabb.min(), m_aabb.max() - m_aabb.min(), m_nodes);
    }

    void forEachDTreeWrapperParallel(std::function<void(DTreeWrapper *)> func) {
        int nDTreeWrappers = static_cast<int>(m_nodes.size());

#pragma omp parallel for
        for (int i = 0; i < nDTreeWrappers; ++i) {
            if (m_nodes[i].isLeaf) {
                func(&m_nodes[i].dTree);
            }
        }
    }

    void record(const Point &p, const Vector &dTreeVoxelSize, DTreeRecord rec,
                EDirectionalFilter directionalFilter,
                EBsdfSamplingFractionLoss bsdfSamplingFractionLoss) {
        float volume = 1;
        for (int i = 0; i < 3; ++i) {
            volume *= dTreeVoxelSize[i];
        }

        rec.statisticalWeight /= volume;
        m_nodes[0].record(p - dTreeVoxelSize * 0.5f,
                          p + dTreeVoxelSize * 0.5f,
                          m_aabb.min(),
                          m_aabb.max() - m_aabb.min(),
                          rec,
                          directionalFilter,
                          bsdfSamplingFractionLoss,
                          m_nodes);
    }

    bool shallSplit(const STreeNode &node, int depth, size_t samplesRequired) {
        return m_nodes.size() < std::numeric_limits<uint32_t>::max() - 1 &&
               node.dTree.statisticalWeightBuilding() > samplesRequired;
    }

    void refine(size_t sTreeThreshold, int maxMB) {
        if (maxMB >= 0) {
            size_t approxMemoryFootprint = 0;
            for (const auto &node : m_nodes) {
                approxMemoryFootprint +=
                    node.dTreeWrapper()->approxMemoryFootprint();
            }

            if (approxMemoryFootprint / 1000000 >= (size_t) maxMB) {
                return;
            }
        }

        struct StackNode {
            size_t index;
            int depth;
        };

        std::stack<StackNode> nodeIndices;
        nodeIndices.push({ 0, 1 });
        while (!nodeIndices.empty()) {
            StackNode sNode = nodeIndices.top();
            nodeIndices.pop();

            // Subdivide if needed and leaf
            if (m_nodes[sNode.index].isLeaf) {
                if (shallSplit(
                        m_nodes[sNode.index], sNode.depth, sTreeThreshold)) {
                    subdivide((int) sNode.index, m_nodes);
                }
            }

            // Add children to stack if we're not
            if (!m_nodes[sNode.index].isLeaf) {
                const STreeNode &node = m_nodes[sNode.index];
                for (int i = 0; i < 2; ++i) {
                    nodeIndices.push({ node.children[i], sNode.depth + 1 });
                }
            }
        }

        // Uncomment once memory becomes an issue.
        // m_nodes.shrink_to_fit();
    }

    const Bounds &aabb() const { return m_aabb; }

private:
    std::vector<STreeNode> m_nodes;
    Bounds m_aabb;
};
}; // namespace lightwave
