#include <lightwave.hpp>

namespace lightwave {

class Sphere : public Shape {
    void populate(SurfaceEvent &surf, const Point &position) const {
        surf.position = Point(Vector(position).normalized());
        surf.uv       = Vector2(
            (std::atan2(-surf.position.z(), surf.position.x()) + Pi) * Inv2Pi,
            safe_acos(surf.position.y()) * InvPi);

        surf.geometryNormal = Vector(surf.position);
        surf.shadingNormal  = surf.geometryNormal;
        surf.tangent        = Vector(1, 0, 0);
        surf.pdf            = Inv4Pi;
    }

public:
    Sphere(const Properties &properties) {}

    bool intersect(const Ray &ray, Intersection &its,
                   Sampler &rng, Context &cont) const override {
        PROFILE("Sphere")
        
        const float b = 2 * ray.direction.dot(Vector(ray.origin));
        const float c = Vector(ray.origin).lengthSquared() - 1;

        float t, tFar;
        if (!solveQuadratic(b, c, t, tFar)) return false;
        if (t < Epsilon) t = tFar;
        if (t < Epsilon || t >= its.t) return false;

        its.t = t;
        cont.curvesIndex = -1;
        populate(its, ray(its.t));
        return true;
    }

    Bounds getBoundingBox() const override {
        return Bounds(Point(-1), Point(+1));
    }

    Point getCentroid() const override { return Point(0); }

    AreaSample sampleArea(Sampler &rng, Context &cont) const override {
        AreaSample sample;
        populate(sample, squareToUniformSphere(rng.next2D()));
        return sample;
    }

    std::string toString() const override { return "Sphere[]"; }
};

} // namespace lightwave

REGISTER_SHAPE(Sphere, "sphere")
