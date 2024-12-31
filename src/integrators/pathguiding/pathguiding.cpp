#include "pathguidingTrees.hpp"
#include <lightwave.hpp>

namespace lightwave {

class PathGuider : public SamplingIntegrator {

    int m_maxDepth;
    bool m_nee;
    bool m_mis;

    float m_bsdfSamplingFraction;
    std::unique_ptr<STree> m_sdTree;
    bool m_isBuilt = false; // SDTree
    int m_budget;
    bool m_useTimeAsBudget = false;
    int m_startingSeconds  = 4;
    int m_sppPerPass       = 4;
    int m_passesRendered   = 0;
    int m_iter             = 0;
    bool m_isFinalIter     = false;
    float m_dTreeThreshold = 0.1f; // subdivide DTree leaf if the power of
                                   // the light is more than this threshold
    int m_sTreeThreshold  = 4000;
    int m_sdTreeMaxMemory = -1; // -1 means unlimited memory usage
    float m_roughnessThreshold; // roughness threshold for using DTree
    std::string m_baseId;
    ESpatialFilter m_spatialFilter;
    EDirectionalFilter m_directionalFilter;
    EBsdfSamplingFractionLoss m_bsdfSamplingFractionLoss;
    EBsdfSamplingFractionLoss m_bsdfSamplingFractionLosssos;
    ENee m_neeType;

    float m_splittingMin;
    float m_splittingMax;
    int m_rrDepth;
    int m_currentRRDepth;
    bool m_debug;

    enum DebugImage {
        None,
        BsdSamplingFraction,
        LREstimate,
        SplittingFactor, // add variables you want to print out here
        All
    } m_debugImage;

    enum ERRSMode {
        ENone,
        EClassic,
        EADRRS,
    };
    ERRSMode m_rrsMethod;
    ERRSMode m_currentRRSmethod = ERRSMode::ENone;

    struct RRSstats {
        int terminations          = 0;
        int splits                = 0;
        int largestSplit          = 0;
        float lowestSplittingVal  = Infinity;
        float highestSplittingVal = 0.0f;
        float totalSplitPaths     = 0;
    };

    RRSstats m_rrsStats;

    void prepareRenderPasses(Context &cont) {
        auto start = std::chrono::high_resolution_clock::now();

        // calculate the amount of samples per pixel in each render pass
        size_t sampleCount         = (size_t) m_budget;
        int nPasses                = std::bit_width(sampleCount) - 1;
        bool result                = true;

        while (result && m_passesRendered < nPasses) {
            int budgetThisIteration =
                m_useTimeAsBudget ? m_startingSeconds << m_iter : 1 << m_iter;

            m_isFinalIter = nPasses == m_passesRendered + 1;

            resetSDTree();

            if (m_neeType == ENee::EAlways ||
                (m_neeType == ENee::EKickstart &&
                 nPasses > m_passesRendered * 2)) {
                m_nee = true;
            } else {
                m_nee = false;
            }
            // logger(EInfo, "NEE: %d", m_nee);

            // train ADRRS with Classic mode
            if (m_debug) {
                if (m_rrsMethod == ERRSMode::EADRRS &&
                    (nPasses > m_passesRendered * 2 || m_passesRendered < 3)) {
                    m_currentRRSmethod = ERRSMode::EClassic;
                    m_currentRRDepth   = 5;
                    logger(EInfo, "Training ADRRS.");
                } else {
                    m_currentRRSmethod = m_rrsMethod;
                    m_currentRRDepth   = m_rrDepth;
                    logger(EInfo, "No longer Training ADRRS.");
                }
                logger(EInfo,
                       "current method: %d",
                       m_currentRRSmethod == ERRSMode::EClassic ? 1
                       : m_currentRRSmethod == ERRSMode::EADRRS ? 2
                                                                : 0);
            }

            m_useTimeAsBudget
                ? performRenderPassSeconds(budgetThisIteration, cont)
                : performRenderPassSpp(budgetThisIteration, cont);

            buildSDTree();

            ++m_iter;
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        logger(EInfo, "Render time: %f", duration.count());
        if (!m_debug) {
            return;
        }
        logger(EInfo,
               "RRS-Statistics: \n Path terminations: %d \n Split paths: %d "
               "\n "
               "Largest split: %d \n Lowest splitting Value: %f \n Highest "
               "Splitting Value: %f \n Average amount of split paths: %f \n",
               m_rrsStats.terminations,
               m_rrsStats.splits,
               m_rrsStats.largestSplit,
               m_rrsStats.lowestSplittingVal,
               m_rrsStats.highestSplittingVal,
               m_rrsStats.totalSplitPaths / m_rrsStats.splits);
    }

    void resetSDTree() {
        logger(EInfo, "Resetting distributions for sampling.");
        /*
        Leaf nodes of the spatial binary tree are subdivided if the number
        of samples they received in the last iteration exceeds c * sqrt(2^k)
        where c is this value and k is the iteration index. The first
        iteration has k==0.
        */
        m_sdTree->refine((size_t) safe_sqrt(pow(2, m_iter) * m_sppPerPass / 4) *
                             m_sTreeThreshold,
                         m_sdTreeMaxMemory);
        m_sdTree->forEachDTreeWrapperParallel([this](DTreeWrapper *dTree) {
            dTree->reset(20, m_dTreeThreshold);
        });
    }

    void buildSDTree() {
        logger(EInfo, "Building distributions for sampling.");

        // build distributions
        m_sdTree->forEachDTreeWrapperParallel(
            [](DTreeWrapper *dTree) { dTree->build(); });
        m_isBuilt = true;

        // everything below is statistics for the logger

        int maxDepth               = 0;
        int minDepth               = std::numeric_limits<int>::max();
        float avgDepth             = 0;
        float maxAvgRadiance       = 0;
        float minAvgRadiance       = std::numeric_limits<float>::max();
        float avgAvgRadiance       = 0;
        size_t maxNodes            = 0;
        size_t minNodes            = std::numeric_limits<size_t>::max();
        float avgNodes             = 0;
        float maxStatisticalWeight = 0;
        float minStatisticalWeight = std::numeric_limits<float>::max();
        float avgStatisticalWeight = 0;

        int nPoints      = 0;
        int nPointsNodes = 0;

        m_sdTree->forEachDTreeWrapperConst([&](const DTreeWrapper *dTree) {
            const int depth = dTree->depth();
            maxDepth        = std::max(maxDepth, depth);
            minDepth        = std::min(minDepth, depth);
            avgDepth += depth;

            const float avgRadiance = dTree->meanRadiance();
            maxAvgRadiance          = max(maxAvgRadiance, avgRadiance);
            minAvgRadiance          = min(minAvgRadiance, avgRadiance);
            avgAvgRadiance += avgRadiance;

            if (dTree->numNodes() > 1) {
                const size_t nodes = dTree->numNodes();
                maxNodes           = std::max(maxNodes, nodes);
                minNodes           = std::min(minNodes, nodes);
                avgNodes += nodes;
                ++nPointsNodes;
            }

            const float statisticalWeight = dTree->statisticalWeight();
            maxStatisticalWeight = max(maxStatisticalWeight, statisticalWeight);
            minStatisticalWeight = min(minStatisticalWeight, statisticalWeight);
            avgStatisticalWeight += statisticalWeight;

            ++nPoints;
        });

        if (nPoints > 0) {
            avgDepth /= nPoints;
            avgAvgRadiance /= nPoints;

            if (nPointsNodes > 0) {
                avgNodes /= nPointsNodes;
            }

            avgStatisticalWeight /= nPoints;
        }

        logger(EInfo,
               "Distribution statistics:\n"
               "  Depth         = [%d, %f, %d]\n"
               "  Mean radiance = [%f, %f, %f]\n"
               "  Node count    = [%zu, %f, %zu]\n"
               "  Stat. weight  = [%f, %f, %f]\n",
               minDepth,
               avgDepth,
               maxDepth,
               minAvgRadiance,
               avgAvgRadiance,
               maxAvgRadiance,
               minNodes,
               avgNodes,
               maxNodes,
               minStatisticalWeight,
               avgStatisticalWeight,
               maxStatisticalWeight);
    }

    void performRenderPassSpp(int samplesPerPixel, Context &cont) {
        logger(
            EInfo, "Render pass with %d samples per pixel.", samplesPerPixel);
        assert_condition((m_image), {
            logger(EError,
                   "<integrator /> needs an <image /> child to render into!");
        });
        m_bsdfSamplingFractionLoss = m_isBuilt
                                         ? m_bsdfSamplingFractionLosssos
                                         : EBsdfSamplingFractionLoss::ENone;

        const Vector2i resolution = m_scene->camera()->resolution();

        const float norm = 1.0f / samplesPerPixel;

        m_image->initialize(resolution);
        if (!m_isFinalIter) {
            m_image->setId(m_baseId + " Iteration: " + std::to_string(m_iter) +
                           " SPP: " + std::to_string(samplesPerPixel));
        } else {
            m_image->setId(m_baseId + " Final Iteration     SPP: " +
                           std::to_string(samplesPerPixel));
        }

        Streaming stream{ *m_image };
        ProgressReporter progress{ resolution.product() };
        for_each_parallel(
            BlockSpiral(resolution, Vector2i(64)), [&](auto block) {
                auto sampler = m_sampler->clone();
                for (auto pixel : block) {
                    Color sum;

                    for (int sample = 0; sample < samplesPerPixel; sample++) {
                        sampler->seed(pixel, sample);
                        auto cameraSample =
                            m_scene->camera()->sample(pixel, *sampler);
                        Color li = Li(cameraSample.ray, *sampler, cont);
                        sum += cameraSample.weight * li;
                    }
                    m_image->get(pixel) = norm * sum;
                }

                progress += block.diagonal().product();
                stream.updateBlock(block);
            });
        progress.finish();

        if (m_isFinalIter) {
            stream.update();
            m_image->save();
        }
        ++m_passesRendered;

        if (!m_debug || m_debugImage == DebugImage::None) {
            return;
        }

        if (m_debugImage != DebugImage::All) {

            Image SingleAOV = Image{ m_image->resolution() };
            SingleAOV.initialize(resolution);
            SingleAOV.setId(m_baseId + " aov " +
                            std::to_string(samplesPerPixel));

            Streaming aovStream{ SingleAOV };
            for (int i = 0; i < resolution.x(); i++) {
                for (int j = 0; j < resolution.y(); j++) {
                    Point2i pixel = Point2i(i, j);
                    auto sampler  = m_sampler->clone();
                    Color sum;

                    sampler->seed(pixel, 0);
                    auto cameraSample =
                        m_scene->camera()->sample(pixel, *sampler);

                    SingleAOV.get(pixel) =
                        aov(cameraSample.ray, *sampler, m_debugImage, cont);
                }
                aovStream.updateBlock(
                    Bounds2i(Point2i(i, 0), Point2i(i + 1, resolution.y())));
            }
        } else {
            Image bsdfSamplingFractionAOV = Image{ m_image->resolution() };
            Image lrEstimateAOV           = Image{ m_image->resolution() };
            Image firstSplitAOV           = Image{ m_image->resolution() };
            bsdfSamplingFractionAOV.initialize(resolution);
            lrEstimateAOV.initialize(resolution);
            firstSplitAOV.initialize(resolution);
            bsdfSamplingFractionAOV.setId(m_baseId +
                                          " bsdfSamplingFraction aov " +
                                          std::to_string(samplesPerPixel));
            lrEstimateAOV.setId(m_baseId + " lrEstimateAOV aov " +
                                std::to_string(samplesPerPixel));
            firstSplitAOV.setId(m_baseId + " firstSplitAOV aov " +
                                std::to_string(samplesPerPixel));
            Streaming aovStream1{ bsdfSamplingFractionAOV };
            Streaming aovStream2{ lrEstimateAOV };
            Streaming aovStream3{ firstSplitAOV };

            for (int i = 0; i < resolution.x(); i++) {
                for (int j = 0; j < resolution.y(); j++) {
                    Point2i pixel = Point2i(i, j);
                    auto sampler  = m_sampler->clone();
                    Color sum;

                    sampler->seed(pixel, 0);
                    auto cameraSample =
                        m_scene->camera()->sample(pixel, *sampler);

                    bsdfSamplingFractionAOV.get(pixel) =
                        aov(cameraSample.ray,
                            *sampler,
                            DebugImage::BsdSamplingFraction,
                            cont);
                    lrEstimateAOV.get(pixel) = aov(cameraSample.ray,
                                                   *sampler,
                                                   DebugImage::LREstimate,
                                                   cont);
                    firstSplitAOV.get(pixel) = aov(cameraSample.ray,
                                                   *sampler,
                                                   DebugImage::SplittingFactor,
                                                   cont);
                }
                aovStream1.updateBlock(
                    Bounds2i(Point2i(i, 0), Point2i(i + 1, resolution.y())));
                aovStream2.updateBlock(
                    Bounds2i(Point2i(i, 0), Point2i(i + 1, resolution.y())));
                aovStream3.updateBlock(
                    Bounds2i(Point2i(i, 0), Point2i(i + 1, resolution.y())));
            }
        }
    }

    void performRenderPassSeconds(int seconds, Context &cont) {
        logger(EInfo, "Render pass for %d sseconds.", seconds);
        assert_condition((m_image), {
            logger(EError,
                   "<integrator /> needs an <image /> child to render into!");
        });
        m_bsdfSamplingFractionLoss = m_isBuilt
                                         ? m_bsdfSamplingFractionLosssos
                                         : EBsdfSamplingFractionLoss::ENone;

        const Vector2i resolution = m_scene->camera()->resolution();

        m_image->initialize(resolution);
        if (!m_isFinalIter) {
            m_image->setId(m_baseId + " Iteration: " + std::to_string(m_iter) +
                           " allocated time: " + std::to_string(seconds));
        } else {
            m_image->setId(m_baseId + " Final Iteration     allocated time: " +
                           std::to_string(seconds));
        }
        int spp = 0;
        Streaming stream{ *m_image };
        auto start = std::chrono::high_resolution_clock::now();
        while (true) {
            spp++;
            stream.normalize(1.0f / spp);
            for_each_parallel(
                BlockSpiral(resolution, Vector2i(64)), [&](auto block) {
                    auto sampler = m_sampler->clone();
                    for (auto pixel : block) {

                        sampler->seed(pixel, spp);
                        auto cameraSample =
                            m_scene->camera()->sample(pixel, *sampler);
                        Color li = Li(cameraSample.ray, *sampler, cont);

                        m_image->get(pixel) += cameraSample.weight * li;
                    }
                    stream.updateBlock(block);
                });

            auto time = std::chrono::high_resolution_clock::now();
            if (std::chrono::duration_cast<std::chrono::seconds>(time - start)
                    .count() >= seconds) {
                break;
            }
        }

        // apply norm at end (not visible in tev)
        const float norm = 1.0f / spp;

        for (int i = 0; i < resolution.x(); i++) {
            for (int j = 0; j < resolution.y(); j++) {
                m_image->get(Point2i(i, j)) *= norm;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        logger(EInfo,
               "Actual render time this iteration: %f\n spp: %d",
               (end - start).count(),
               spp);

        if (m_isFinalIter) {
            stream.update();
            m_image->save();
        }
        ++m_passesRendered;
    }

    Color aov(const Ray &ray, Sampler &rng, DebugImage aovType, Context &cont) {
        Intersection its = m_scene->intersect(ray, rng, cont);
        Vector dTreeVoxelSize;
        DTreeWrapper *dTree =
            m_sdTree->dTreeWrapper(its.position, dTreeVoxelSize);
        switch (aovType) {
        case DebugImage::BsdSamplingFraction: {
            return Color(dTree->bsdfSamplingFraction());
        }
        case DebugImage::LREstimate: {
            return Color(dTree->getLrEstimate());
        }
        case DebugImage::SplittingFactor: {
            if (!its.instance || !its.instance->bsdf()) {
                return Color(0);
            }
            return Color(evaluateRR(dTree,
                                    its.instance->bsdf()->albedo(its.uv, cont),
                                    Color(1),
                                    1));
        }
        default:
            break;
        }
        return Color(0);
    }

    float weightWindow(float splittingFactor,
                       float weightWindowSize = 5) const {
        const float dminus = 2 / (1 + weightWindowSize);
        const float dplus  = dminus * weightWindowSize;

        if (splittingFactor < dminus) {
            /// russian roulette
            return splittingFactor / dminus;
        } else if (splittingFactor > dplus) {
            /// splitting
            return splittingFactor / dplus;
        } else {
            /// within weight window
            return 1;
        }
    }

    float evaluateRR(const DTreeWrapper *dTree, const Color &albedo,
                     const Color &throughput, int depth) const {
        if (depth < m_currentRRDepth) {
            /// do not perform RR(S) at this depth.
            return 1.0f;
        }

        switch (m_currentRRSmethod) {
        case ERRSMode::EADRRS: {
            /// "Adjoint-driven Russian Roulette and Splitting"
            Color LiEstimate = dTree->getLrEstimate();
            if (LiEstimate.r() > 0 || LiEstimate.g() > 0 ||
                LiEstimate.b() > 0) {
                return clamp(weightWindow((throughput * LiEstimate).mean()),
                             m_splittingMin,
                             m_splittingMax);
            } else {
                return clamp(1.0f, m_splittingMin, m_splittingMax);
            }
        }
        case ERRSMode::EClassic:
            /// Classic RR based on throughput weight
            if (albedo.r() == 0.0f && albedo.g() == 0.0f && albedo.b() == 0.0f)
                /// avoid bias for materials that might report their
                /// reflectance incorrectly
                return 0.1f;
            return clamp((throughput * albedo).mean(), 0.0f, 0.95f);
        case ERRSMode::ENone:
            return clamp(1.0f, m_splittingMin, m_splittingMax);
        }
        return 0;
    }

    float miWeight(float a, float b) const {
        if (std::isinf(a))
            return 1;
        return a / (a + b);
    }

    Color nextEventEstimation(const Intersection &its, const float time,
                              Sampler &rng, const DTreeWrapper *dTree,
                              float bsdfSamplingFraction,
                              const Vector dTreeVoxelSize,
                              Context &cont) const {
        PROFILE("NEE")

        auto lightSample = m_scene->sampleLight(rng);
        if (!lightSample)
            return Color(0);

        auto directSample = lightSample.light->sampleDirect(its.position, rng);
        if (m_scene->intersect(Ray(its.position, directSample.wi, time, 0),
                               directSample.distance,
                               rng,
                               cont))
            // occluded
            return Color(0);

        auto B = its.evaluateBsdf(directSample.wi, cont);

        float misWeight = 1.0f;
        if (m_mis && lightSample.light->canBeIntersected() && B) {
            float woPdf;
            if (!dTree || !m_isBuilt ||
                its.instance->bsdf()->getRoughness(its.uv, cont) <=
                    m_roughnessThreshold) {
                woPdf = B.pdf;
            } else if (!std::isfinite(B.pdf)) {
                woPdf = 0.0f;
            } else {
                float dTreePdf = dTree->pdf(directSample.wi);
                woPdf          = bsdfSamplingFraction * B.pdf +
                        (1 - bsdfSamplingFraction) * dTreePdf;
            }
            misWeight =
                miWeight(lightSample.probability * directSample.pdf, woPdf);
        }

        return B.value * misWeight * directSample.weight /
               lightSample.probability;
    }

    // samples according to m_normal bsdf or dtree
    BsdfSample sampleMaterial(Sampler &rng, const DTreeWrapper *dTree,
                              const Intersection &its, float &woPdf,
                              float &dTreePdf, const float bsdfSamplingFraction,
                              Context &cont) {
        if (!its.instance || !its.instance->bsdf()) {
            return BsdfSample::invalid();
        }

        BsdfSample result;

        if (!dTree || !m_isBuilt ||
            its.instance->bsdf()->getRoughness(its.uv, cont) <= m_roughnessThreshold) {
            result   = its.sampleBsdf(rng, cont);
            woPdf    = result.pdf;
            dTreePdf = 0;
            return result;
        }

        float sample = rng.next();
        Color colorValue;
        float bsdfPdf;

        if (sample < bsdfSamplingFraction) {

            BsdfSample bsdfSamp = its.sampleBsdf(rng, cont);
            colorValue          = bsdfSamp.weight * bsdfSamp.pdf;

            bsdfPdf = bsdfSamp.pdf;

            if (bsdfSamp.isInvalid()) {
                woPdf = dTreePdf = 0;
                return bsdfSamp;
            }

            if (std::isinf(bsdfSamp.pdf)) {
                dTreePdf = 0;
                woPdf    = Infinity;
                bsdfSamp.weight /= bsdfSamplingFraction;

                return bsdfSamp;
            }
            result.wi = bsdfSamp.wi;
        } else {
            Vector wi         = dTree->sample(&rng);
            BsdfEval bsdfEval = its.evaluateBsdf(wi, cont);

            bsdfPdf    = bsdfEval.pdf;
            colorValue = bsdfEval.value;
            result.wi  = wi;
        }

        // finalize pdf

        dTreePdf = dTree->pdf(result.wi);
        woPdf    = bsdfSamplingFraction * bsdfPdf +
                (1 - bsdfSamplingFraction) * dTreePdf;

        result.pdf    = bsdfPdf;
        result.weight = colorValue / woPdf;
        return result;
    }

    // write the power of the lightpath at each Vertex in the correct DTree
    void splatLightPoint(DTreeWrapper *dTree, const float woPdf,
                         const float dTreePdf, const float bsdfPdf,
                         const Color bsdfValue, const Color radiance,
                         const Ray ray, const Vector dTreeVoxelSize,
                         Sampler &rng) {

        if (!(woPdf > 0) || !radiance.isValid() || !bsdfValue.isValid()) {
            return;
        }

        Color product = radiance * bsdfValue;
        float statisticalWeight =
            m_neeType == ENee::EKickstart && m_nee ? 0.5f : 1.0f;
        bool isDelta = bsdfPdf == Infinity;

        DTreeRecord rec = { ray.direction,     radiance.mean(),
                            product.mean(),    woPdf,
                            bsdfPdf,           dTreePdf,
                            statisticalWeight, isDelta };

        switch (m_spatialFilter) {
        case ESpatialFilter::ENearest:
            dTree->record(rec, m_directionalFilter, m_bsdfSamplingFractionLoss);
            break;
        case ESpatialFilter::EStochasticBox: {
            DTreeWrapper *splatDTree = dTree;
            Vector offset            = dTreeVoxelSize;
            offset *= rng.next3D() - Vector(0.5f);
            Point origin = m_sdTree->aabb().clip(ray.origin + offset);
            splatDTree   = m_sdTree->dTreeWrapper(origin);
            if (splatDTree) {
                splatDTree->record(
                    rec, m_directionalFilter, m_bsdfSamplingFractionLoss);
            }

            break;
        }
        case ESpatialFilter::EBox:
            m_sdTree->record(ray.origin,
                             dTreeVoxelSize,
                             rec,
                             m_directionalFilter,
                             m_bsdfSamplingFractionLoss);
            break;
        }
    }

    Color recursiveLi(const Ray &ray, Sampler &rng, int depth, float prevWoPdf,
                      Color throughput, Context &cont) {
        auto L   = Color(0);
        auto its = m_scene->intersect(ray, rng, cont);
        if (auto E = its.evaluateEmission(cont)) {
            const auto misWeight =
                std::isinf(prevWoPdf) || !m_nee || (E.pdf == 0) ? 1
                : m_mis ? miWeight(prevWoPdf, E.pdf)
                        : 0;
            L += misWeight * E.value;
        }

        if (!its || (ray.depth >= m_maxDepth)) {
            return L;
        }

        Vector dTreeVoxelSize;
        DTreeWrapper *dTree =
            m_sdTree->dTreeWrapper(its.position, dTreeVoxelSize);

        float bsdfSamplingFraction = m_bsdfSamplingFraction;

        if (!its.instance || !its.instance->bsdf()) {
            return L;
        }
        const Color albedo = its.instance->bsdf()->albedo(its.uv, cont);
        Color res          = Color(0);
        float splittingFactor =
            evaluateRR(dTree, albedo, throughput, ray.depth);
        const int numSamples = int(splittingFactor + rng.next());

        //  prepare stats for rr debugging
        if (m_debug) {
            if (numSamples == 0) {
                m_rrsStats.terminations++;
            } else if (numSamples > 1) {
                m_rrsStats.splits++;
                m_rrsStats.totalSplitPaths += numSamples;
            }
            m_rrsStats.largestSplit = max(numSamples, m_rrsStats.largestSplit);
            m_rrsStats.lowestSplittingVal =
                min(m_rrsStats.lowestSplittingVal, splittingFactor);
            m_rrsStats.highestSplittingVal =
                max(m_rrsStats.highestSplittingVal, splittingFactor);
        }

        for (int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex) {

            if (dTree && m_bsdfSamplingFractionLoss !=
                             EBsdfSamplingFractionLoss::ENone) {
                bsdfSamplingFraction = dTree->bsdfSamplingFraction();
            }

            if (m_nee) {
                // MISweight already applied in nextEventEstimation
                L += nextEventEstimation(its,
                                         ray.time,
                                         rng,
                                         dTree,
                                         bsdfSamplingFraction,
                                         dTreeVoxelSize,
                                         cont);
            }

            // assuming all bsdfs are smooth
            float woPdf, dTreePdf;
            auto M = sampleMaterial(
                rng, dTree, its, woPdf, dTreePdf, bsdfSamplingFraction, cont);

            if (!M)
                return L;

            Color woWeight = M.weight;

            Ray newRay = Ray(its.position, M.wi, ray.time, ray.depth + 1);
            Color incomingLight =
                recursiveLi(newRay,
                            rng,
                            depth + 1,
                            woPdf,
                            throughput * woWeight / splittingFactor,
                            cont);

            if (dTree && 1 / woPdf > 0) {
                // FIXME: hier auch splitting factor rein?

                splatLightPoint(dTree,
                                woPdf,
                                dTreePdf,
                                M.pdf,
                                M.weight / woPdf,
                                incomingLight,
                                newRay,
                                dTreeVoxelSize,
                                rng);
            }

            res = (L + woWeight * incomingLight);

            if (!res.isValid()) {
                logger(EInfo,
                       "Res: %s\n L: %s\nwoWeight: %s\nincomingLight: %s",
                       res,
                       L,
                       woWeight,
                       incomingLight);
            }
        }
        // splat power (only necessary when ADRRS will be used)
        if (m_rrsMethod == ERRSMode::EADRRS && numSamples > 0) {
            dTree->splatLrEstimate(res, 1.f * numSamples);
        }
        return res / splittingFactor;
    }

public:
    PathGuider(const Properties &properties) : SamplingIntegrator(properties) {
        m_maxDepth = properties.get<int>("depth", 2);
        m_mis      = properties.get<bool>("mis", false);
        m_bsdfSamplingFraction =
            properties.get<float>("bsdfSamplingFraction", 0.5f);
        m_budget          = properties.get<int>("budget", 1024);
        m_useTimeAsBudget = properties.get<bool>("timeAsBudget", false);
        m_startingSeconds = properties.get<int>("startingSeconds", 4);
        m_roughnessThreshold = m_bsdfSamplingFraction == 1.0f ? 1.0f :
            properties.get<float>("roughnessThreshold", 0.2f);
        m_neeType =
            properties.getEnum<ENee>("neeType",
                                     ENee::EAlways,
                                     { { "kickstart", ENee::EKickstart },
                                       { "always", ENee::EAlways },
                                       { "never", ENee::ENever } });
        m_spatialFilter = properties.getEnum<ESpatialFilter>(
            "spatialFilter",
            ESpatialFilter::ENearest,
            { { "nearest", ESpatialFilter::ENearest },
              { "box", ESpatialFilter::EBox },
              { "stochasticBox", ESpatialFilter::EStochasticBox } });
        m_directionalFilter = properties.getEnum<EDirectionalFilter>(
            "directionalFilter",
            EDirectionalFilter::ENearest,
            { { "nearest", EDirectionalFilter::ENearest },
              { "box", EDirectionalFilter::EBox } });
        m_bsdfSamplingFractionLosssos =
            properties.getEnum<EBsdfSamplingFractionLoss>(
                "bsdfSamplingFractionLoss",
                EBsdfSamplingFractionLoss::ENone,
                { { "none", EBsdfSamplingFractionLoss::ENone },
                  { "kl", EBsdfSamplingFractionLoss::EKL },
                  { "variance", EBsdfSamplingFractionLoss::EVariance } });
        m_rrsMethod =
            properties.getEnum<ERRSMode>("rrsMethod",
                                         ERRSMode::ENone,
                                         { { "none", ERRSMode::ENone },
                                           { "classic", ERRSMode::EClassic },
                                           { "adrrs", ERRSMode::EADRRS } });
        m_splittingMin = properties.get<float>("splittingMin", 0.05);
        m_splittingMax = properties.get<float>("splittingMax", 20.0f);
        m_rrDepth      = properties.get<int>("rrDepth", 1);
        m_debug        = properties.get<bool>("debug", false);
        m_debugImage   = properties.getEnum<DebugImage>(
            "debugAOV",
            DebugImage::None,
            { { "bsdfSamplingFraction", DebugImage::BsdSamplingFraction },
                { "lrEstimate", DebugImage::LREstimate },
                { "splittingFactor", DebugImage::SplittingFactor },
                { "all", DebugImage::All },
                { "none", DebugImage::None } });

        if (m_rrsMethod == ERRSMode::EADRRS && m_rrDepth != 2) {
            logger(EWarn, "ADRRS should ideally be used with rrDepth 2");
        }
        if (m_rrsMethod == ERRSMode::EClassic && m_rrDepth != 5) {
            logger(EWarn, "Classic RR should ideally be used with rrDepth 5");
        }
    }

    void execute() override {
        // TODO: check if this is correct
        Context cont;
        m_sdTree = std::unique_ptr<STree>(new STree(scene()->getBoundingBox()));
        m_baseId = m_image->id();
        prepareRenderPasses(cont);
    }

    Color Li(const Ray &cameraRay, Sampler &rng, Context &cont) override {
        PROFILE("Li")

        return recursiveLi(cameraRay, rng, 1, Infinity, Color(1), cont);
    }

    std::string toString() const override {
        return tfm::format(
            "PathTracer[\n"
            "  sampler = %s,\n"
            "  image = %s,\n"
            "]",
            indent(m_sampler),
            indent(m_image));
    }
};

} // namespace lightwave

REGISTER_INTEGRATOR(PathGuider, "pathguider")
