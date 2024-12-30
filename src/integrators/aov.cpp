#include <lightwave.hpp>

namespace lightwave {

class AovIntegrator : public SamplingIntegrator {
    enum Variable {
        AovNormals,
        AovDistance,
        AovAlbedo,
        AovBvh,
        AovUv,
    } m_variable;
    float m_scale;

public:
    AovIntegrator(const Properties &properties)
        : SamplingIntegrator(properties) {
        // clang-format off
        m_variable = properties.getEnum<Variable>("variable", {
            { "normals",  AovNormals  },
            { "distance", AovDistance },
            { "albedo",   AovAlbedo   },
            { "bvh",      AovBvh      },
            { "uv",       AovUv       },
        });
        m_scale = properties.get<float>("scale", 1.f);
        // clang-format on
    }

    Color Li(const Ray &ray, Sampler &rng, Context &cont) override {
        Intersection its = m_scene->intersect(ray, rng, cont);

        switch (m_variable) {
        case AovNormals:
            return its ? (Color(its.shadingNormal) + Color(1)) / 2
                       : Color(0.5f);
        case AovDistance:
            return its ? Color(its.t) / m_scale : Color(Infinity);
        case AovAlbedo:
            if (its)
                if (auto bsdf = its.instance->bsdf())
                    return bsdf->albedo(its.uv, cont);
            return Color(0);
        case AovBvh:
            return Color(its.stats.bvhCounter / m_scale,
                         its.stats.primCounter / m_scale,
                         0);
        case AovUv:
            return its ? Color(its.uv.x(), its.uv.y(), 0) : Color(0);
        default:
            return Color(0);
        }
    }

    std::string toString() const override {
        return tfm::format(
            "AovIntegrator[\n"
            "  sampler = %s,\n"
            "  image = %s,\n"
            "]",
            indent(m_sampler),
            indent(m_image));
    }
};

} // namespace lightwave

REGISTER_INTEGRATOR(AovIntegrator, "aov")
