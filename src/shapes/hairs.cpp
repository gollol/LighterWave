#include <bob.hpp>
#include <lightwave.hpp>

#include "accel.hpp"

namespace lightwave {

// See
// https://research.nvidia.com/sites/default/files/pubs/2018-08_Phantom-Ray-Hair-Intersector//Phantom-HPG%202018.pdf
/// @brief We expect the curve's control points to be passed as Bezier.
class Hairs : public AccelerationStructure {

protected:
    // Almost identically taken from the paper.
    struct RayConeIntersection { // ray.origin = Point{0}, ray.direction =
                                 // Vector{0,0,1}
        Point c0;                // curve(t) in RCC
        Vector cd;               // tangent(t) in RCC
        float s;                 // the intersection's t paramter
        float dt;                // unscaled delta(t) from the paper
        float dp;
        float dc;
        float sp;

        RayConeIntersection(Point pointOnCurve, Vector tangentAtPoint) {
            c0 = pointOnCurve;
            cd = tangentAtPoint;
        }

        // Returns true for real intersections and false for phantom
        // intersections.
        inline bool intersect(const float r, const float dr) {
            const float r2  = r * r;
            const float drr = r * dr;

            float ddd = sqr(cd.x()) + sqr(cd.y());
            dp        = sqr(c0.x()) + sqr(c0.y());

            const float cdd = c0.x() * cd.x() +
                              c0.y() * cd.y(); // Has to be (c(t) -
                                               // o).dot(c'(t)) according to sp.
            const float cxd  = c0.x() * cd.y() - c0.y() * cd.x();
            const float cdz2 = sqr(cd.z());

            float c = ddd;
            ddd += cdz2;
            float b = cd.z() * (drr - cdd);
            float a = 2 * drr * cdd + cxd * cxd - ddd * r2 + dp * cdz2;

            const float det = b * b - a * c;
            s               = (b - safe_sqrt(det)) / c;
            dt              = (s * cd.z() - cdd) / ddd;
            dc              = s * s + dp;
            sp              = cdd / cd.z();
            dp += sp * sp;

            s += c0.z();
            sp += c0.z();
            return det > 0;
        }
    };

    // Parametric version of a cylinder 'tightly' fitting the hair,
    // used for early exit checks
    struct BoundingCylinder {
        Vector de;
        Point oe;
        float maxDist;
    };

    std::vector<BoundingCylinder> m_boundingCylinders; // Precomputed
    std::vector<Bounds> m_boundingBoxes; // Will be cleared after creating the
                                         // acceleration structure
    std::vector<Point> m_controlPoints;  // Contains ALL control points in order

    // These three vectors are only used to keep track of indizes when
    // reordering. They slow the program down by quite a bit (~ 20% for some
    // scenes) and are pretty large but other solutions are WAY too ugly.
    std::vector<int> m_hairIndizes;
    std::vector<int> m_curveIndizes;
    std::vector<int> m_splitIndizes;

    std::filesystem::path m_originalPath;
    float m_rootRadius;
    float m_tipRadius;

    // The number of control points per strand
    int32_t m_numControlPoints;
    // The number of segments per strand
    int32_t m_numCurves;
    // The number of strands
    int32_t m_hairCount;

    // Number of sampled points when computing bounding boxes and cylinders
    static constexpr int32_t m_numCurveSamples = 20;
    // Max number of iterations when looking for an intersection
    static constexpr int32_t m_phantomIter = 30;
    // Number of sub-segments each segment is split into (paper uses 8)
    static constexpr int32_t m_curveSplits = 16;

    static constexpr float m_invCurveSplits = 1.f / m_curveSplits;
    static constexpr float m_margin         = 1e-5;
    static constexpr float m_deltaMargin    = 5e-5;

protected:
    int numberOfPrimitives() const override {
        return m_hairCount * m_numCurves * m_curveSplits;
    }

    Bounds getBoundingBox(int primitiveIndex) const override {
        return m_boundingBoxes[primitiveIndex];
    }

    /// @brief Gets the centroid of the bounding box as there is no clear
    /// centroid of a curve
    Point getCentroid(int primitiveIndex) const override {
        return m_boundingBoxes[primitiveIndex].center();
    }

    /// @brief t in range [0, 1] and curveIndex = 0 on first bezier curve, 1 on
    /// second, ...
    inline float localRadius(const float t, const int curveIndex) const {
        const float hairT = (t + curveIndex) / m_numCurves;
        return interpolateLinear(hairT, m_rootRadius, m_tipRadius);
        // return interpolateWave(hairT, 9, m_rootRadius, m_tipRadius);
    }

    /// @brief t in range [0, 1] and curveIndex = 0 on first bezier curve, 1 on
    /// second, ...
    inline float localSlant(const float t, const int curveIndex) const {
        const float hairT = (t + curveIndex) / m_numCurves;
        return interpolateLinearPrime(hairT, m_rootRadius, m_tipRadius) /
               m_numCurves;
        // return interpolateWavePrime(hairT, 9, m_rootRadius, m_tipRadius) /
        // m_numCurves;
    }

    /// @brief Contains the 'high-level' workflow when looking for an
    /// intersection
    bool intersect(const int primitiveIndex, const Ray &ray, Intersection &its,
                   Sampler &rng, Context &cont) const override {
        const int splitIndex = m_splitIndizes[primitiveIndex];
        const int curveIndex = m_curveIndizes[primitiveIndex];
        const int hairIndex  = m_hairIndizes[primitiveIndex];
        const int cpIndex    = hairIndex * m_numControlPoints + 3 * curveIndex;

        const float tstart = splitIndex * m_invCurveSplits;
        const float tend   = (splitIndex + 1) * m_invCurveSplits;

        if (!intersectsBoundingCylinder(tstart, primitiveIndex, ray))
            return false;

        const Point w0               = m_controlPoints[cpIndex];
        const Point w1               = m_controlPoints[cpIndex + 1];
        const Point w2               = m_controlPoints[cpIndex + 2];
        const Point w3               = m_controlPoints[cpIndex + 3];
        const std::array<Point, 4> w = { w0, w1, w2, w3 };

        SimpleTransform rccToLocal = getRccTransform(ray, w3);

        // // If deltaT(start) < 0 and deltaT(end) > 0, ignore interval
        if (delT(tstart, curveIndex, w, rccToLocal) < 0 &&
            delT(tend, curveIndex, w, rccToLocal) > 0)
            return false;

        const Point start = interpolateBezierCurve(tstart, w);
        const Point end   = interpolateBezierCurve(tend, w);

        // Start iterations at 'better' end and afterwards, try the 'worse' one.
        const float t1 = ray.direction.dot(end - start) > 0 ? tstart : tend;
        const float t2 = ray.direction.dot(end - start) > 0 ? tend : tstart;
        if (intersectionIteration(
                t1, curveIndex, hairIndex, ray, its, cont, rccToLocal, w))
            return true;
        if (intersectionIteration(
                t2, curveIndex, hairIndex, ray, its, cont, rccToLocal, w))
            return true;

        return false;
    }

    /// @brief Returns false if the ray misses the bounding cylinder of the
    /// curve.
    inline bool intersectsBoundingCylinder(const float tstart,
                                           const int primitiveIndex,
                                           const Ray &ray) const {
        const BoundingCylinder bCylinder = m_boundingCylinders[primitiveIndex];
        const Vector n                   = bCylinder.de.cross(ray.direction);

        const float distRayCylindAxis = abs(n.dot(ray.origin - bCylinder.oe));
        return distRayCylindAxis < bCylinder.maxDist;
    }

    /// @brief Returns the transform going from ray-centric coordinates to local
    /// space
    inline SimpleTransform getRccTransform(const Ray &ray,
                                           const Point &w3) const {
        const Vector d = ray.direction;
        Vector q       = Vector(w3).cross(d);

        // If w3 and d are (nearly) collinear we choose an arbitrary q
        q = (q.length() < m_margin) ? Frame(d).tangent : q.normalized();
        const Vector c = q.cross(d);

        Matrix3x3 changeOfBasis; // This is an orthonormal matrix!
        changeOfBasis.setColumn(0, q);
        changeOfBasis.setColumn(1, c);
        changeOfBasis.setColumn(2, d);

        StaticTransform rccToLoc(changeOfBasis, true);
        rccToLoc.translate(Vector(ray.origin));

        // Using rccToLoc.inverse(), we transform into a coordinate
        // system where ray.direction = {0,0,1} and ray.origin = {0,0,0}

        return rccToLoc.interpolate(0);
    }

    /// @brief Contains the entire root-finding logic (i.e.
    /// intersection-finding) and returns true iff an intersection is found.
    /// Starting at t = tstart, we use the RayConeIntersection struct from the
    /// paper to update t until we're either close enough or our iterations run
    /// dry.
    inline bool intersectionIteration(const float tstart, const int curveIndex,
                                      const int hairIndex, const Ray &ray,
                                      Intersection &its, Context &cont,
                                      const SimpleTransform &rccToLocal,
                                      const std::array<Point, 4> &cp) const {
        float t          = tstart;
        float prevDeltaT = 0, prevT;
        for (int k = 0; k < m_phantomIter; k++) {
            const Point rccStart =
                rccToLocal.inverse(interpolateBezierCurve(t, cp));
            const Vector rccTangent =
                rccToLocal.inverse((Vector) tangentToBezierCurve(t, cp));
            RayConeIntersection coneInter(rccStart, rccTangent);

            const bool realRoot = coneInter.intersect(
                localRadius(t, curveIndex), localSlant(t, curveIndex));
            const float deltaT = clamp(coneInter.dt, -0.5f, 0.5f);

            prevT = t;
            if (k == 0) { // We want to do atleast 2 iterations
                prevDeltaT = deltaT;
                continue;
            }

            // When we have a 0-crossing for delta(t), "it's more stable to do
            // regula falsi", otherwise we just add deltaT
            if (deltaT * prevDeltaT >= 0)
                t += deltaT;
            else {
                if (k % 4 == 0)
                    t = (t + prevT) / 2;
                else
                    t = (deltaT * prevT - prevDeltaT * t) /
                        (deltaT - prevDeltaT);
            }

            // Buttend handling
            const float buffer = min(0.1, 2 * m_invCurveSplits);
            const float lower  = max(0.f, tstart - buffer);
            const float upper  = min(1.f, tstart + m_invCurveSplits + buffer);
            if (t < lower || t > upper) {
                // If we neither hit the root nor the tip, then the 'buttend'
                // lies between 2 segments and is thus occluded.
                if (!(curveIndex == 0 && t < 0.f) &&
                    !(curveIndex == m_numCurves - 1 && t > 1.f))
                    return false;

                const bool tbool =
                    (t >= 1.f) ? true
                               : false; // Safely cast t into a bool. (True if t
                                        // > 1 and False if t < 0)
                t = clamp(t, 0.f, 1.f);
                const Vector tangent =
                    Vector(tangentToBezierCurve(t, cp)).normalized();
                const Point endPoint = tbool ? cp[3] : cp[0];

                const float denom = tangent.dot(ray.direction);
                if (abs(denom) < m_margin)
                    return false;

                const float tbar = (endPoint - ray.origin).dot(tangent) / denom;
                if (tbar < Epsilon || tbar > its.t)
                    return false;

                const Point hitPoint = ray(tbar);
                if (abs((hitPoint - endPoint).length()) >=
                    localRadius(t, curveIndex))
                    return false;

                const Vector gNormal = (tbool ? 1.f : -1.f) * tangent;

                cont.curvesIndex   = hairIndex;
                cont.tangentNormal = gNormal;
                cont.hairThickness = localRadius(t, curveIndex);
                cont.hairIntercept = (t + curveIndex) / m_numCurves;

                its.t  = tbar;
                its.uv = Point2{ 1.f * curveIndex / (m_numCurves - 1), 0.f };
                its.position       = hitPoint;
                its.geometryNormal = gNormal;
                its.shadingNormal  = its.geometryNormal;
                its.tangent        = hitPoint - endPoint;
                its.pdf            = 0;

                assert_finite(its.t,
                              { logger(EError, "offending shape: %s", this); });
                return true;
            }

            // Stopping the iteration when we are precise enough
            if (abs(deltaT) < m_deltaMargin) {
                const float s = coneInter.s;
                if (!realRoot || s < Epsilon || s >= its.t)
                    return false;

                if (abs(deltaT - prevDeltaT) > m_margin)
                    t = (deltaT * prevT - prevDeltaT * t) /
                        (deltaT - prevDeltaT);

                const float hairT   = (t + curveIndex) / m_numCurves;
                const Point onCurve = interpolateBezierCurve(t, cp);
                const Point onCone  = ray(s);
                const Vector n      = onCone - onCurve;
                const Vector gNormal =
                    (n.length() > 1e-8 ? n.normalized() : -ray.direction);
                const Vector tangentCurve =
                    Vector(tangentToBezierCurve(t, cp)).normalized();
                const Vector tanBitan = tangentCurve.cross(-ray.direction);
                Vector sNormal        = gNormal;

                // Try this one if the geometric normal doesn't look right!
                // sNormal = (gNormal +
                // 0.875f*gNormal.dot(ray.direction)*ray.direction).normalized();

                // https://github.com/blender/blender/blob/ddbc34829feb0ab46b2723d87e6602d76117b9c7/intern/cycles/kernel/geom/curve_intersect.h
                // Generates the same normals as Curve Info -> Tangent Normals,
                // no idea wtf that is though.
                cont.tangentNormal = -tangentCurve.cross(tanBitan);
                cont.curvesIndex   = hairIndex;
                cont.hairThickness = localRadius(t, curveIndex);
                cont.hairIntercept = hairT;

                its.t  = s;
                its.uv = Point2{ hairT, 0 }; // Don't know how to consistently
                                             // calculate v :(
                its.pdf            = 0;
                its.position       = onCone;
                its.geometryNormal = gNormal; // Is -ray.direction for Curve
                                              // type 'Rounded Ribbons'
                its.shadingNormal = sNormal;
                its.tangent =
                    (tangentCurve - tangentCurve.dot(gNormal) * gNormal)
                        .normalized(); // Projection of the tangent to the curve
                                       // onto the cone

                assert_finite(its.t,
                              { logger(EError, "offending shape: %s", this); });
                return true;
            }

            prevDeltaT = deltaT;
        }
        return false;
    }

    /// @brief Calculates the delta(t) from the paper if that's all we need
    inline float delT(const float t, int curveIndex,
                      const std::array<Point, 4> &w,
                      const SimpleTransform &rccToLocal) const {
        const Point start = rccToLocal.inverse(interpolateBezierCurve(t, w));
        const Vector tangent =
            rccToLocal.inverse((Vector) tangentToBezierCurve(t, w));

        RayConeIntersection coneInter(start, tangent);
        coneInter.intersect(localRadius(t, curveIndex),
                            localSlant(t, curveIndex));
        return coneInter.dt;
    }

protected:
    /// @brief Computes the bounding box and cylinder for each hair segment
    inline void populateBoundingObjects() {
        for (int hairIndex = 0; hairIndex < m_hairCount; hairIndex++) {
            int cpIndex      = hairIndex * m_numControlPoints;
            int segmentIndex = hairIndex * m_numCurves;

            for (int i = 0; i < m_numControlPoints - 1; i += 3) {
                const Point w0               = m_controlPoints[cpIndex];
                const Point w1               = m_controlPoints[cpIndex + 1];
                const Point w2               = m_controlPoints[cpIndex + 2];
                const Point w3               = m_controlPoints[cpIndex + 3];
                const std::array<Point, 4> w = { w0, w1, w2, w3 };

                // We split each curve segment into parts to find all delta(T)
                // solutions and get a better bounding box fit
                Point start, end = interpolateBezierCurve(0, w);
                for (int numSplit = 0; numSplit < m_curveSplits; numSplit++) {
                    start = end;
                    end   = interpolateBezierCurve(
                        (numSplit + 1.f) * m_invCurveSplits, w);

                    Bounds aabb = Bounds::empty();
                    const Point oe =
                        0.5f *
                        (start + interpolateBezierCurve(
                                     (numSplit + 0.5f) * m_invCurveSplits, w));
                    const Vector de = (end - start).normalized();
                    float maxDist   = 0;

                    // We sample points on the curve to approximate the max
                    // distance from the curve to the line 'oe + lambda * de'
                    // Instead of using the control points, we use the sampled
                    // points for a tighter bounding box
                    for (int j = 0; j < m_numCurveSamples; j++) {
                        float t = (m_numCurveSamples > 1)
                                      ? (1.f * j) / (m_numCurveSamples - 1)
                                      : 0;
                        t       = (numSplit + t) * m_invCurveSplits;

                        // Building bounding box
                        const float radius = localRadius(t, i / 3);
                        const Point sample = interpolateBezierCurve(t, w);
                        aabb.extend(sample + Vector(radius));
                        aabb.extend(sample - Vector(radius));

                        // Sampling bounding cylinder's maxDist
                        const Vector toLine = oe - interpolateBezierCurve(t, w);
                        const float distInterpolToLine =
                            (toLine - toLine.dot(de) * de).length();
                        if (distInterpolToLine > maxDist)
                            maxDist = distInterpolToLine;
                    }
                    maxDist += localRadius(numSplit * m_invCurveSplits, i / 3);

                    m_boundingCylinders[m_curveSplits * segmentIndex +
                                        numSplit] = { .de      = de,
                                                      .oe      = oe,
                                                      .maxDist = maxDist };
                    m_boundingBoxes[m_curveSplits * segmentIndex + numSplit] =
                        aabb;
                }
                cpIndex += 3;
                segmentIndex++;
            }
        }
    }

    void loadCurves(const std::filesystem::path &path) {
        bob::Fur fur(path);
        const int furSize = fur.getKeyPoints().size();
        assert_condition(
            furSize >= 4,
            logger(EInfo, "Need atleast 4 control points, got %d", furSize););

        m_numControlPoints = fur.getKeyPointCount();
        m_controlPoints.reserve(furSize);

        for (auto p : fur.getKeyPoints())
            m_controlPoints.push_back({ p.data() });

        const int cpSize = m_controlPoints.size();
        assert_condition(furSize == cpSize,
                         logger(EInfo,
                                "Hair was loaded but not correctly transfered, "
                                "got a difference of %d",
                                furSize - cpSize););
        assert_condition(m_numControlPoints % 3 == 1,
                         logger(EInfo,
                                "We need 3n+1 many control points for cubic "
                                "Bezier splines, got %d",
                                m_numControlPoints););
        assert_condition(cpSize % m_numControlPoints == 0,
                         logger(EInfo,
                                "Not all hair strands have the same amout of "
                                "control points!"););

        m_numCurves = (m_numControlPoints - 1) / 3;
        m_hairCount = cpSize / m_numControlPoints;
    }

    void reorderPrimitives(const std::vector<int> &newOrder) override {
        std::vector<BoundingCylinder> updatedCylinders;
        std::vector<Bounds> updatedBoundingBoxes;
        std::vector<int> updatedHairIndizies;
        std::vector<int> updatedCurveIndizies;
        std::vector<int> updatedSplitIndizies;
        for (int index : newOrder) {
            updatedCylinders.push_back(m_boundingCylinders[index]);
            updatedBoundingBoxes.push_back(m_boundingBoxes[index]);
            updatedHairIndizies.push_back(m_hairIndizes[index]);
            updatedCurveIndizies.push_back(m_curveIndizes[index]);
            updatedSplitIndizies.push_back(m_splitIndizes[index]);
        }
        m_boundingCylinders = updatedCylinders;
        m_boundingBoxes     = updatedBoundingBoxes;
        m_hairIndizes       = updatedHairIndizies;
        m_curveIndizes      = updatedCurveIndizies;
        m_splitIndizes      = updatedSplitIndizies;
    }

public:
    Hairs(const Properties &properties) {
        m_originalPath = properties.get<std::filesystem::path>("filename");
        m_rootRadius   = properties.get<float>("rootRadius");
        m_tipRadius    = properties.get<float>("tipRadius");

        loadCurves(m_originalPath);
        logger(EInfo,
               "loaded bob with %d strands, %d control points",
               m_hairCount,
               m_controlPoints.size());

        m_boundingBoxes =
            std::vector<Bounds>(m_hairCount * m_numCurves * m_curveSplits);
        m_boundingCylinders = std::vector<BoundingCylinder>(
            m_hairCount * m_numCurves * m_curveSplits);

        m_hairIndizes =
            std::vector<int>(m_hairCount * m_numCurves * m_curveSplits);
        m_curveIndizes =
            std::vector<int>(m_hairCount * m_numCurves * m_curveSplits);
        m_splitIndizes =
            std::vector<int>(m_hairCount * m_numCurves * m_curveSplits);

        for (int i = 0; i < m_hairCount * m_numCurves * m_curveSplits; i++) {
            m_hairIndizes[i] = i / (m_numCurves * m_curveSplits);
            m_curveIndizes[i] =
                (i % (m_numCurves * m_curveSplits)) / m_curveSplits;
            m_splitIndizes[i] = i % m_curveSplits;
        }

        populateBoundingObjects();

        buildAccelerationStructure();
        m_boundingBoxes.clear();
    }

    bool intersect(const Ray &ray, Intersection &its, Sampler &rng,
                   Context &cont) const override {
        PROFILE("Hair")
        return AccelerationStructure::intersect(ray, its, rng, cont);
    }

    AreaSample sampleArea(Sampler &rng,
                          Context &cont) const override{ NOT_IMPLEMENTED }

    std::string toString() const override {
        return tfm::format(
            "Hairs[\n"
            "  control points = %s,\n"
            "  strands = %s,\n"
            "  filename = \"%s\"\n"
            "]",
            indent(m_controlPoints.size()),
            indent(m_hairCount),
            indent(m_originalPath));
    }
};

} // namespace lightwave

REGISTER_SHAPE(Hairs, "fur")
