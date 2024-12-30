#include <lightwave.hpp>

namespace lightwave {

class Mandelbulb : public Shape {
    static constexpr float Power      = 3;
    static constexpr float SdfEpsilon = 1e-4;
    static constexpr int MaxSteps     = 200;

    float smoothmin(float a, float b, float strength) const {
        return -log(exp(-a * strength) + exp(-b * strength)) / strength;
    }

    float sdf(const Point &pos, int &steps) const {
        // float a = (pos - Point(-0.5,0,0)).length() - 1;
        // float b = (pos - Point(+0.5,0,0)).length() - 1;
        // float k = 20;
        // return smoothmin(a, b, 10);
        // return min(a, b);

        Vector z{ pos };
        float dr = 1;
        float r  = 0;
        for (steps = 0; steps < 20; steps++) {
            r = z.length();
            if (r > 4)
                break;

            // convert to polar coordinates
            float theta = acos(z.z() / r);
            float phi   = atan2(z.y(), z.x());
            dr          = pow(r, Power - 1) * Power * dr + 1;

            // scale and rotate the point
            float zr = pow(r, Power);
            theta    = theta * Power;
            phi      = phi * Power;

            // convert back to cartesian coordinates
            z = zr * Vector(sin(theta) * cos(phi), sin(theta) * sin(phi),
                            cos(theta)) +
                Vector(pos);
        }

        return 0.5f * std::log(r) * r / dr;
    }

    Vector normal(const Point &pos) const {
        int steps;
        return (Vector(sdf(pos, steps)) -
                Vector(sdf(pos + Vector(SdfEpsilon, 0, 0), steps),
                       sdf(pos + Vector(0, SdfEpsilon, 0), steps),
                       sdf(pos + Vector(0, 0, SdfEpsilon), steps)))
            .normalized();
    }

public:
    Mandelbulb(const Properties &properties) {}

    bool intersect(const Ray &ray, Intersection &its,
                   Sampler &rng, Context &cont) const override {
        int steps;

        float t = Epsilon;
        for (int i = 0; i < MaxSteps; i++) {
            float dist = sdf(ray(t), steps);
            if (dist >= its.t || dist > 1e+3)
                break;
            if (dist < SdfEpsilon) {
                its.t = t;

                its.position = ray(t);
                its.uv.x()   = its.position.x();
                its.uv.y()   = its.position.y();
                its.shadingNormal = normal(its.position);
                its.geometryNormal = normal(its.position);
                its.tangent = Vector(1, 0, 0);

                cont.curvesIndex = -1;
                return true;
            }
            t += dist;
        }

        return false;
    }

    Bounds getBoundingBox() const override {
        return Bounds(Point(-10), Point(+10));
    }

    Point getCentroid() const override { return Point(0); }

    std::string toString() const override { return "Mandelbulb[]"; }
};

} // namespace lightwave

REGISTER_SHAPE(Mandelbulb, "mandelbulb")
