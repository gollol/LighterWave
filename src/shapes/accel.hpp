#pragma once

#include <lightwave/core.hpp>
#include <lightwave/math.hpp>
#include <lightwave/shape.hpp>

#include <bitset>
#include <numeric>

namespace lightwave {

/**
 * @brief Parent class for shapes that combine many individual shapes (e.g.,
 * triangle meshes), and hence benefit from building an acceleration structure
 * over their children.
 *
 * To use this class, you will need to implement the following methods:
 * - numberOfPrimitives()           -- report the number of individual children
 * that the shape has
 * - intersect(primitiveIndex, ...) -- intersect a single child (identified by
 * the given index) for the given ray
 * - getBoundingBox(primitiveIndex) -- return the bounding box of a single child
 * (used for building the BVH)
 * - getCentroid(primitiveIndex)    -- return the centroid of a single child
 * (used for building the BVH)
 *
 * @example For a simple example of how to use this class, look at @ref
 * shapes/group.cpp
 * @see Group
 * @see TriangleMesh
 */
class AccelerationStructure : public Shape {
    /// @brief The datatype used to index BVH nodes and the primitive index
    /// remapping.
    typedef int32_t NodeIndex;

    /// @brief The number of bits used per compressed point component.
    static constexpr uint32_t cBits = 12;
    /// @brief The biggest representable value for a point component.
    static constexpr uint32_t c = static_cast<uint32_t>((1 << cBits) - 1);
    /// @brief The inverse of c.
    static constexpr float cInv = 1.0f / c;

    /// @brief The extends of the uncompressed root node aka the diagonal vector
    /// of the BVH.
    Vector m_extends;
    /// @brief The inverse of the diagonal vector of the root node (1 /
    /// diagonal).
    Vector m_invExtends;
    /// @brief The lower bound of the uncompressed root node.
    Vector m_lowerBound;

    /// @brief A point with each component compressed to a c-bit integer.
    struct CompressedPoint {

        enum Position { Lower, Upper };
        unsigned int x : cBits;
        unsigned int y : cBits;
        unsigned int z : cBits;

        CompressedPoint() = default;
        CompressedPoint(const Vector &p, Position pos) {
            if (pos == Lower) {
                x = (uint32_t) floor(clamp(p.x(), 0.0f, 1.0f) * c);
                y = (uint32_t) floor(clamp(p.y(), 0.0f, 1.0f) * c);
                z = (uint32_t) floor(clamp(p.z(), 0.0f, 1.0f) * c);
            } else {
                x = (uint32_t) ceil(clamp(p.x(), 0.0f, 1.0f) * c);
                y = (uint32_t) ceil(clamp(p.y(), 0.0f, 1.0f) * c);
                z = (uint32_t) ceil(clamp(p.z(), 0.0f, 1.0f) * c);
            }
        }

        /// @brief Decompresses the point
        Vector decompress() const {
            return { static_cast<float>(x) * cInv,
                     static_cast<float>(y) * cInv,
                     static_cast<float>(z) * cInv };
        }

        /// @brief Returns the compressed point with integer components
        operator Vector() const {
            return { static_cast<float>(x),
                     static_cast<float>(y),
                     static_cast<float>(z) };
        }
    };

    /// @brief A compressed version of Bounds using compressed points as min and
    /// max positions.
    struct BoundsStorage {
        CompressedPoint low, high;
        BoundsStorage() = default;
        BoundsStorage(CompressedPoint low, CompressedPoint high)
            : low(low), high(high) {}
        BoundsStorage(const Bounds &b, const Vector &invExtends,
                      const Point &lowerBound) {
            PROFILE("Compress bounds");
            Vector lowerDist = (b.min() - lowerBound) * invExtends;
            Vector upperDist = (b.max() - lowerBound) * invExtends;

            low  = CompressedPoint(lowerDist, CompressedPoint::Lower);
            high = CompressedPoint(upperDist, CompressedPoint::Upper);
        }

        /// @brief  Returns the compressed bounding box of the root node
        static BoundsStorage root() {
            return { { Vector(0), CompressedPoint::Lower },
                     { Vector(1), CompressedPoint::Upper } };
        }

        /// @brief  Returns the empty bounding box
        static BoundsStorage empty() {
            return { { Vector(0), CompressedPoint::Lower },
                     { Vector(0), CompressedPoint::Lower } };
        }

        /// @brief returns the raw, compressed data of the bounding volume
        operator Bounds() const { return { low, high }; }

        /// @brief returns the bounding volume decompressed to Bounds
        Bounds decompress(const Vector &extends,
                          const Vector &lowerBound) const {
            PROFILE("Decompress bounds");
            return { low.decompress() * extends + lowerBound,
                     high.decompress() * extends + lowerBound };
        }
    };

    /// @brief A node in our binary BVH tree.
    struct Node {
        /// @brief The axis aligned bounding box of this node.
        BoundsStorage aabb;
        /**
         * @brief Either the index of the left child node in the bvh (for
         * internal nodes), or the first primitive in m_primitiveIndices (for
         * leaf nodes). Primitive indices are stored as negative numbers,
         * which is also used to distinguish internal nodes from leaf nodes.
         */
        NodeIndex left;
        /// @brief The index of the right child node in the bvh,
        /// or the number of primitives for leaf nodes.
        NodeIndex right;

        static Node invalid() { return { BoundsStorage::empty(), 0, 0 }; }

        /// @brief Whether this node is invalid.
        bool isInvalid() const { return left == 0 && right == 0; }

        /// @brief Whether this BVH node is a leaf node.
        bool isLeaf() const { return left <= 0; }

        /// @brief For internal nodes: The index of the left child node in
        /// the bvh.
        NodeIndex leftChildIndex() const { return left; }
        /// @brief For internal nodes: The index of the right child node in
        /// the bvh.
        NodeIndex rightChildIndex() const { return right; }

        /// @brief For leaf nodes: The first index in m_primitiveIndices.
        NodeIndex firstPrimitiveIndex() const { return -left; }
        /// @brief For leaf nodes: The last index in m_primitiveIndices (still
        /// included).
        NodeIndex lastPrimitiveIndex() const { return -left + right - 1; }

        /// @brief The number of primitives in this node.
        NodeIndex primitiveCount() const { return right; }

        /// @brief Set the number of primitives in this node.
        void setPrimitiveCount(int count) { right = count; }
    };

    struct BuilderNode : Node {
        Bounds dc_aabb;
        /// @brief Marks the node as the root node of a cluster.
        bool isClusterRoot = false;
        /// @brief Marks the node as the left child of its parent.
        bool isLeftChild = false;
        /// @brief Marks the node as the right child of its parent.
        bool isRightChild = false;

        /// @brief The size of the subtree rooted at this node.
        int subtreeSize = 0;

        /// @brief The index of the parent node in the bvh.
        NodeIndex parentIndex = 0;

        /// @brief Create a node that holds tow child-clusters after a merge
        static BuilderNode ghost() {
            BuilderNode n;
            n.isLeftChild  = true;
            n.isRightChild = true;
            return n;
        }

        /// @brief Whether this node is a ghost node.
        bool isGhost() const { return isLeftChild && isRightChild; }
        /// @brief Whether this node is the root of the whole BVH-tree.
        bool isTreeRoot() const { return !(isLeftChild || isRightChild); }
        /// @brief Marks this node as the root of the whole BVH-tree.
        void markTreeRoot() {
            isLeftChild  = false;
            isRightChild = false;
        }
    };

    /// @brief A ray in BoundsStorage-space with inverted direction (1 / dir)
    /// and other useful data for intersecting the BVH in compressed space
    struct TransformedRay {
        Point origin;
        Vector invDirection;
        float length;
        float t;

        TransformedRay(const Ray &r, Intersection &its,
                       const AccelerationStructure &as) {
            origin = (r.origin - Point(as.m_lowerBound)) * as.m_invExtends *
                     static_cast<float>(c);
            Vector dir = r.direction * as.m_invExtends * static_cast<float>(c);

            length       = dir.length();
            invDirection = Vector(1.0f) / (dir / length);
            t            = its.t * length;
        }
    };

    /// @brief A list of all BVH nodes.
    std::vector<Node> m_nodes;

    /// @brief Returns the root BVH node.
    const Node &rootNode() const {
        // by convention, this is always the first element of m_nodes
        return m_nodes.front();
    }

    /**
     * @brief Intersects a BVH node, recursing into children (for internal
     * nodes), or intersecting all primitives (for leaf nodes).
     */
    bool intersectNode(const Node &node, const Ray &original, Intersection &its,
                       Sampler &rng, TransformedRay &tr, Context &cont) const {
        // update the statistic tracking how many BVH nodes have been tested
        // for intersection
        its.stats.bvhCounter++;

        bool wasIntersected = false;
        if (node.isLeaf()) {
            for (NodeIndex i = 0; i < node.primitiveCount(); i++) {
                // test the child for intersection
                wasIntersected |= intersect(
                    node.firstPrimitiveIndex() + i, original, its, rng, cont);
            }
            // update the statistic tracking how many children have been
            // tested for intersection
            its.stats.primCounter += node.primitiveCount();
            tr.t = its.t * tr.length;
        } else { // internal node
            // test which bounding box is intersected first by the ray.
            // this allows us to traverse the children in the order they are
            // intersected in, which can help prune a lot of unnecessary
            // intersection tests.
            const Node &left  = m_nodes[node.leftChildIndex()];
            const Node &right = m_nodes[node.rightChildIndex()];

            const auto leftT  = intersectAABB(left.aabb, tr);
            const auto rightT = intersectAABB(right.aabb, tr);
            if (leftT < rightT) { // left child is hit first; test left
                                  // child
                // first, then right child
                if (leftT < tr.t)
                    wasIntersected |=
                        intersectNode(left, original, its, rng, tr, cont);
                if (rightT < tr.t)
                    wasIntersected |=
                        intersectNode(right, original, its, rng, tr, cont);
            } else { // right child is hit first; test right child first,
                     // then
                // left child
                if (rightT < tr.t)
                    wasIntersected |=
                        intersectNode(right, original, its, rng, tr, cont);
                if (leftT < tr.t)
                    wasIntersected |=
                        intersectNode(left, original, its, rng, tr, cont);
            }
        }
        return wasIntersected;
    }

    /// @brief Performs a slab test to intersect a bounding box with a ray,
    /// returning Infinity in case the ray misses.
    float intersectAABB(const BoundsStorage &sBounds,
                        const TransformedRay &tr) const {
        Bounds bounds = sBounds;
        const auto t1 = (bounds.min() - tr.origin) * tr.invDirection;
        // intersect all axes at once with the maximum slabs of the bounding
        // box
        const auto t2 = (bounds.max() - tr.origin) * tr.invDirection;

        // the elementwiseMin picks the near slab for each axis, of which we
        // then take the maximum
        const auto tNear = elementwiseMin(t1, t2).maxComponent();
        // the elementwiseMax picks the far slab for each axis, of which we
        // then take the minimum
        const auto tFar = elementwiseMax(t1, t2).minComponent();

        if (tFar < tNear)
            return Infinity; // the ray does not intersect the bounding box
        if (tFar < Epsilon)
            return Infinity; // the bounding box lies behind the ray origin

        return tNear; // return the first intersection with the bounding box
        // (may also be negative!)
    }

protected:
    struct BVHBuilder {
        typedef BuilderNode BNode;
        typedef std::vector<BNode> BTree;

        /**
         * @brief Mapping from internal @c NodeIndex to @c primitiveIndex as
         * used by all interface methods. For efficient storage, we assume
         * that children of BVH leaf nodes have contiguous indices, which
         * would require re-ordering the primitives. For simplicity, we
         * instead perform this re-ordering on a list of indices (which
         * starts of as
         * @code 0, 1, 2, ..., primitiveCount - 1 @endcode ), which allows
         * us to translate from re-ordered (contiguous) indices to the
         * indices the user of this class expects.
         */
        std::vector<int> m_primitiveIndices;
        std::vector<int> m_newOrder;

        AccelerationStructure &as;
        std::vector<BNode> m_heavyNodes;

        BVHBuilder(AccelerationStructure &as) : as(as) {};

        /// @brief Builds and optimizes the bvh for the underlying shape.
        void build() {
            // fill primitive indices with 0 to primitiveCount - 1
            m_primitiveIndices.resize(as.numberOfPrimitives());
            std::iota(m_primitiveIndices.begin(), m_primitiveIndices.end(), 0);

            // create root node
            auto &root = m_heavyNodes.emplace_back();
            root.left  = 0;
            root.setPrimitiveCount(as.numberOfPrimitives());
            root.markTreeRoot();

            // build the bvh-tree
            computeRootAABB(root);
            subdivide(0);

            // cache optimize the bvh-tree
            as.m_nodes.reserve(m_heavyNodes.size());
            m_newOrder.reserve(m_primitiveIndices.size());
            optimize(m_heavyNodes.front());

            as.reorderPrimitives(m_newOrder);
        }

    private:
        /// @brief Computes the axis aligned bounding box for the root node
        void computeRootAABB(BNode &node) {
            node.dc_aabb = Bounds::empty();
            for (NodeIndex i = 0; i < node.primitiveCount(); i++) {
                const Bounds childAABB = as.getBoundingBox(
                    m_primitiveIndices[node.firstPrimitiveIndex() + i]);
                node.dc_aabb.extend(childAABB);
            }

            as.m_lowerBound = Vector(node.dc_aabb.min());
            as.m_extends    = node.dc_aabb.diagonal();

            // safe version to compute inverse extents
            for (size_t i = 0; i < 3; i++) {
                as.m_invExtends[i] =
                    as.m_extends[i] < Epsilon ? 1.0f : 1.0f / as.m_extends[i];
            }

            node.aabb = BoundsStorage::root();
        }

        /// @brief Computes the axis aligned bounding box for a leaf BVH
        /// node
        void computeAABB(BNode &node) {
            node.dc_aabb = Bounds::empty();
            for (NodeIndex i = 0; i < node.primitiveCount(); i++) {
                const Bounds childAABB = as.getBoundingBox(
                    m_primitiveIndices[node.firstPrimitiveIndex() + i]);
                node.dc_aabb.extend(childAABB);
            }
            node.aabb =
                BoundsStorage(node.dc_aabb, as.m_invExtends, as.m_lowerBound);
            node.dc_aabb = node.aabb.decompress(as.m_extends, as.m_lowerBound);
        }

        /// @brief Computes the surface area of a bounding box.
        float surfaceArea(const Bounds &bounds) const {
            const auto size = bounds.diagonal();
            return 2 * (size.x() * size.y() + size.x() * size.z() +
                        size.y() * size.z());
        }

        /// @brief Computes the overlap of two bounding boxes.
        float overlap(const Bounds &a, const Bounds &b) const {
            Bounds overlap = a.clip(b);
            if (overlap.isEmpty())
                return 0.f;
            return surfaceArea(overlap);
        }
        /// @brief Computes the combined surface area of the two children of a
        /// node
        float combinedChildArea(const BNode &parent, BTree &bvh) const {
            if (parent.isLeaf())
                return 0.0f;

            const Bounds &left  = bvh[parent.left].dc_aabb;
            const Bounds &right = bvh[parent.right].dc_aabb;
            return surfaceArea(left) + surfaceArea(right) -
                   overlap(left, right);
        }

        /**
         * For a given node, computes split axis and split position that
         * minimize the surface area heuristic.
         * @param node The BVH node to compute the split for.
         * @param out bestSplitAxis The optimal split axis, or -1 if no
         * useful split exists
         * @param out bestSplitPosition The optimal split position,
         * undefined if no useful split exists
         */
        void binning(const BNode &node, int &bestSplitAxis,
                     float &bestSplitPosition) {
            static constexpr size_t BinCount = 16;
            struct Bin {
                Bounds aabb;
                NodeIndex primCount{ 0 };
                float rightCost{ 0.f };
            };
            struct Split {
                float cost     = Infinity;
                int axis       = -1;
                float position = 0;
            };

            Split bestSplit;

            // compute SAH cost of doing no split and use as baseline
            const float traversalCost = 1.f;
            bestSplit.cost = (node.primitiveCount() - traversalCost) *
                             surfaceArea(node.dc_aabb);

            const NodeIndex firstPrim = node.firstPrimitiveIndex();
            const NodeIndex lastPrim  = node.lastPrimitiveIndex();

            // for (int test = 0; test < 1; test++) {
            // const int splitAxis =
            // node.aabb.diagonal().maxComponentIndex();
            for (int splitAxis = 0; splitAxis < 3; splitAxis++) {
                Split split;
                split.axis = splitAxis;

                // compute range of centroid
                float centroidMin = +Infinity;
                float centroidMax = -Infinity;
                for (NodeIndex i = firstPrim; i <= lastPrim; i++) {
                    const float centroid =
                        as.getCentroid(m_primitiveIndices[i])[splitAxis];
                    centroidMin = std::min(centroidMin, centroid);
                    centroidMax = std::max(centroidMax, centroid);
                }

                if (centroidMin == centroidMax) {
                    continue;
                }

                const float binSize = (centroidMax - centroidMin) / BinCount;
                const float inverseBinSize =
                    BinCount / (centroidMax - centroidMin);

                // compute bins
                std::array<Bin, BinCount> bins;
                for (NodeIndex i = firstPrim; i <= lastPrim; i++) {
                    const float centroid =
                        as.getCentroid(m_primitiveIndices[i])[split.axis];
                    int binId = int((centroid - centroidMin) * inverseBinSize);
                    binId     = min(int(BinCount - 1), max(binId, 0));

                    auto &bin = bins[binId];
                    bin.aabb.extend(as.getBoundingBox(m_primitiveIndices[i]));
                    bin.primCount++;
                }

                // Sweep bins to compute SAH
                Bounds sweepBBox;
                NodeIndex sweepCount = 0;
                for (int i = BinCount - 1; i >= 0; --i) {
                    sweepCount += bins[i].primCount;
                    sweepBBox.extend(bins[i].aabb);
                    bins[i].rightCost = surfaceArea(sweepBBox) * sweepCount;
                }

                sweepBBox  = Bounds{};
                sweepCount = 0;
                for (size_t i = 0; i < BinCount - 1; ++i) {
                    sweepCount += bins[i].primCount;
                    sweepBBox.extend(bins[i].aabb);
                    const float leftCost  = surfaceArea(sweepBBox) * sweepCount;
                    const float totalCost = leftCost + bins[i + 1].rightCost;
                    if (totalCost < split.cost) {
                        split.cost     = totalCost;
                        split.position = centroidMin + (int(i) + 1) * binSize;
                    }
                }

                if (split.cost < bestSplit.cost) {
                    bestSplit = split;
                }
            }

            bestSplitAxis     = bestSplit.axis;
            bestSplitPosition = bestSplit.position;
        }

        /// @brief Attempts to subdivide a given BVH node.
        void subdivide(NodeIndex i) {
            BNode &parent = m_heavyNodes[i];
            // only subdivide if enough children are available.
            if (parent.primitiveCount() <= 2) {
                return;
            }

            static constexpr bool UseSAH = true;
            // set to true when implementing binning

            int splitAxis = -1;
            float splitPosition;
            if (UseSAH) {
                // pick split axis and position using binned SAH
                binning(parent, splitAxis, splitPosition);
            } else {
                // split in the middle of the longest axis
                splitAxis     = parent.dc_aabb.diagonal().maxComponentIndex();
                splitPosition = parent.dc_aabb.center()[splitAxis];
            }

            if (splitAxis == -1) {
                // a split axis of -1 indicates that no useful split exists
                return;
            }

            // the point at which to split (note that primitives must be
            // re-ordered so that all children of the left node will have a
            // smaller index than firstRightIndex, and nodes on the right
            // will have an index larger or equal to firstRightIndex)
            NodeIndex firstRightIndex = parent.firstPrimitiveIndex();
            NodeIndex lastLeftIndex   = parent.lastPrimitiveIndex();

            // partition algorithm (you might remember this from quicksort)
            while (firstRightIndex <= lastLeftIndex) {
                if (as.getCentroid(
                        m_primitiveIndices[firstRightIndex])[splitAxis] <
                    splitPosition) {
                    firstRightIndex++;
                } else {
                    std::swap(m_primitiveIndices[firstRightIndex],
                              m_primitiveIndices[lastLeftIndex--]);
                }
            }

            const NodeIndex firstLeftIndex = parent.firstPrimitiveIndex();
            const NodeIndex leftCount      = firstRightIndex - firstLeftIndex;
            const NodeIndex rightCount = parent.primitiveCount() - leftCount;

            if (leftCount == 0 || rightCount == 0) {
                // if either child gets no primitives, we abort subdividing
                return;
            }

            const NodeIndex leftChildIndex = m_heavyNodes.size();
            m_heavyNodes[i].left           = leftChildIndex;
            m_heavyNodes.emplace_back();
            m_heavyNodes[leftChildIndex].left = -firstLeftIndex;
            m_heavyNodes[leftChildIndex].setPrimitiveCount(leftCount);
            m_heavyNodes[leftChildIndex].isLeftChild = true;
            m_heavyNodes[leftChildIndex].parentIndex = i;

            // first, process the left child node (and all of its children)
            computeAABB(m_heavyNodes[leftChildIndex]);
            subdivide(leftChildIndex);

            const NodeIndex rightChildIndex = m_heavyNodes.size();
            m_heavyNodes[i].right           = rightChildIndex;
            m_heavyNodes.emplace_back();
            m_heavyNodes[rightChildIndex].left = -firstRightIndex;
            m_heavyNodes[rightChildIndex].setPrimitiveCount(rightCount);
            m_heavyNodes[rightChildIndex].isRightChild = true;
            m_heavyNodes[rightChildIndex].parentIndex  = i;

            // then, process the right child node (and all of its children)
            computeAABB(m_heavyNodes[rightChildIndex]);
            subdivide(rightChildIndex);
            m_heavyNodes[i].subtreeSize = m_heavyNodes.size() - i;
        }

        /// @brief Recursively optimizes the bvh in regard to cache locality
        void optimize(BNode &root) {
            if (root.subtreeSize <= 3) {
                storeSubtree(root);
                return;
            }
            std::vector<NodeIndex> childRoots = computeRootCluster(root);

            updateRootClusterSize(root);
            optimize(root);

            while (childRoots.size() > 0) {
                BNode c = mergeChildren(childRoots);
                optimize(c);
            }
        }

        /// @brief Puts nodes of the current cluster in m_nodes
        void storeSubtree(BNode &root) {
            if (root.isGhost()) {
                BNode &left = m_heavyNodes[root.left];
                if (!left.isClusterRoot)
                    storeNodes(left);

                BNode &right = m_heavyNodes[root.right];
                if (!right.isClusterRoot)
                    storeNodes(right);
            } else
                storeNodes(root);
        }

        void storeNodes(BNode &root) {
            if (!root.isTreeRoot()) {
                Node &parent = as.m_nodes[root.parentIndex];
                if (root.isLeftChild)
                    parent.left = as.m_nodes.size();
                else
                    parent.right = as.m_nodes.size();
            }

            int parentIndex = as.m_nodes.size();

            if (root.isLeaf()) {
                int firstPrimitiveIndex = root.firstPrimitiveIndex();
                root.left               = -m_newOrder.size();
                for (int i = 0; i < root.primitiveCount(); i++) {
                    m_newOrder.push_back(
                        m_primitiveIndices[firstPrimitiveIndex + i]);
                }
                as.m_nodes.emplace_back(root);
                return;
            }

            as.m_nodes.emplace_back(root);

            BNode &left      = m_heavyNodes[root.left];
            left.parentIndex = parentIndex;
            if (!left.isClusterRoot)
                storeNodes(left);

            BNode &right      = m_heavyNodes[root.right];
            right.parentIndex = parentIndex;
            if (!right.isClusterRoot)
                storeNodes(right);
        }

        /// @brief Computes which node in the bvh should be part of the root
        /// cluster
        /// @return The children of the root cluster
        std::vector<NodeIndex> computeRootCluster(BNode &root) {
            int clusterSize  = ceil(sqrt(root.subtreeSize + 1.0f) - 1.0f);
            root.subtreeSize = clusterSize;

            std::vector<NodeIndex> children;
            children.reserve(clusterSize);

            if (!m_heavyNodes[root.left].isClusterRoot)
                children.emplace_back(root.left);
            if (!m_heavyNodes[root.right].isClusterRoot)
                children.emplace_back(root.right);

            for (int i = root.isGhost() ? 0 : 1; i < clusterSize; i++) {
                float maxArea        = 0.0f;
                size_t maxChildIndex = 0;

                // find highest probability child
                for (size_t j = 0; j < children.size(); j++) {
                    BNode &child = m_heavyNodes[children[j]];
                    float area   = surfaceArea(child.dc_aabb);
                    if (area > maxArea) {
                        maxArea       = area;
                        maxChildIndex = j;
                    }
                }

                BNode &newNode = m_heavyNodes[children[maxChildIndex]];
                children.erase(children.begin() + maxChildIndex);

                if (newNode.isLeaf())
                    continue;

                if (!m_heavyNodes[newNode.left].isClusterRoot)
                    children.emplace_back(newNode.left);
                if (!m_heavyNodes[newNode.right].isClusterRoot)
                    children.emplace_back(newNode.right);
            }

            // mark leafs of new root cluster
            for (NodeIndex c : children)
                m_heavyNodes[c].isClusterRoot = true;

            return children;
        }

        /// @brief Updates the subtree size of the roots in the root cluster
        void updateRootClusterSize(BNode &root) {
            if (root.isLeaf()) {
                root.subtreeSize = 1;
                return;
            }
            BNode &left  = m_heavyNodes[root.left];
            BNode &right = m_heavyNodes[root.right];
            if (!left.isClusterRoot)
                updateRootClusterSize(left);
            if (!right.isClusterRoot)
                updateRootClusterSize(right);

            root.subtreeSize = !left.isClusterRoot * left.subtreeSize +
                               !right.isClusterRoot * right.subtreeSize + 1;
        }

        /// @brief Merges two child-clusters into one in case of a strong
        /// overlap. If the conditions for a merge are not met, a child is
        /// returned.
        BNode mergeChildren(std::vector<NodeIndex> &childRoots) {
            NodeIndex childIndex = childRoots.back();
            BNode &child         = m_heavyNodes[childIndex];
            childRoots.pop_back();

            float childrenArea = combinedChildArea(child, m_heavyNodes);

            NodeIndex bestCousinIndex = -1;
            float maxOverlapRoots     = 0.0f;

            for (size_t j = 0; j < childRoots.size(); j++) {
                BNode &cousin      = m_heavyNodes[childRoots[j]];
                float overlapRoots = overlap(cousin.dc_aabb, child.dc_aabb);

                if (overlapRoots > childrenArea &&
                    overlapRoots > maxOverlapRoots) {
                    float siblingArea = combinedChildArea(cousin, m_heavyNodes);

                    if (overlapRoots > siblingArea) {
                        bestCousinIndex = j;
                        maxOverlapRoots = overlapRoots;
                    }
                }
            }

            if (bestCousinIndex < 0)
                return child;

            NodeIndex cousinIndex = childRoots[bestCousinIndex];
            BNode &cousin         = m_heavyNodes[cousinIndex];
            childRoots.erase(childRoots.begin() + bestCousinIndex);
            child.isClusterRoot  = false;
            cousin.isClusterRoot = false;

            BNode parent       = BNode::ghost();
            parent.left        = cousinIndex;
            parent.right       = childIndex;
            parent.subtreeSize = cousin.subtreeSize + child.subtreeSize;

            return parent;
        }
    };

    /// @brief Returns the number of children (individual shapes) that are
    /// part of this acceleration structure.
    virtual int numberOfPrimitives() const = 0;
    /// @brief Intersect a single child (identified by the index) with the
    /// given ray.
    virtual bool intersect(int primitiveIndex, const Ray &ray,
                           Intersection &its, Sampler &rng,
                           Context &cont) const = 0;
    /// @brief Returns the axis aligned bounding box of the given child.
    virtual Bounds getBoundingBox(int primitiveIndex) const = 0;
    /// @brief Returns the centroid of the given child.
    virtual Point getCentroid(int primitiveIndex) const = 0;
    /// @brief Reorders the internal primitive vector based in the input.
    virtual void reorderPrimitives(const std::vector<int> &newOrder) = 0;

    /// @brief Builds the acceleration structure.
    void buildAccelerationStructure() {
        PROFILE("Build BVH");
        Timer buildTimer;

        BVHBuilder builder(*this);
        builder.build();

        logger(EInfo,
               "built BVH with %ld nodes for %ld primitives in %.1f ms",
               m_nodes.size(),
               numberOfPrimitives(),
               buildTimer.getElapsedTime() * 1000);
    }

public:
    bool intersect(const Ray &ray, Intersection &its, Sampler &rng,
                   Context &cont) const override {
        PROFILE("Accel");
        // TODO: add this back in again (in a way that works)
        // if (m_primitiveIndices.empty())
        //     return false; // exit early if no children exist

        TransformedRay tr(ray, its, *this);

        if (intersectAABB(rootNode().aabb, tr) < tr.t) // test root
                                                       // bounding
                                                       // box for
                                                       // potential
                                                       // hit
            return intersectNode(rootNode(), ray, its, rng, tr, cont);
        return false;
    }

    Bounds getBoundingBox() const override {
        return rootNode().aabb.decompress(m_extends, m_lowerBound);
    }

    Point getCentroid() const override {
        return rootNode().aabb.decompress(m_extends, m_lowerBound).center();
    }
};

} // namespace lightwave
