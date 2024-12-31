#include <lightwave.hpp>

namespace lightwave {


class Spherical : public Lens {

public:
    Spherical(const Properties &properties): Lens(properties) {}

    bool intersect(const Ray &ray, const float zPosition, float &t, Vector &n) const override {

        const float zCenter = zPosition + m_curvatureRadius; 
        const float radiusSquared = sqr(m_curvatureRadius);

        // Calculate quadratic coefficients
        Point oc = ray.origin - Vector(0, 0, zCenter);
        Vector dir = ray.direction;
        float b = 2 * (dir.x() * oc.x() + dir.y() * oc.y() + dir.z() * oc.z());
        float c = sqr(oc.x()) + sqr(oc.y()) + sqr(oc.z()) - radiusSquared;

        // Solve quadratic equation for t values
        float t0, t1;
        if (!solveQuadratic(b, c, t0, t1)) return false;

        // Select intersection t based on ray direction and element curvature
        bool useCloserT = (dir.z() > 0) ^ (m_curvatureRadius < 0);
        t = useCloserT ? std::min(t0, t1) : std::max(t0, t1);
        if (t < 0) return false;

        // Compute surface normal of element at intersection point
        Point pIts = ray(t);
        n = (pIts - Point(0, 0, zCenter)).normalized();
        // Flip normal to point towards incoming ray direction
        if (n.dot(ray.direction) > 0) n = -n;
    
        return true;
    }

    std::string toString() const override {
        return tfm::format(
            "Spherical[\n"
            "  thickness = %d,\n"
            "  eta = %d,\n"
            "  apertureRadius = %d,\n"
            "  curvatureRadius = %d,\n"
            "]",
            m_thickness,
            m_eta,
            m_apertureRadius,
            m_curvatureRadius
        );
    }

};

}

REGISTER_LENS(Spherical, "spherical")