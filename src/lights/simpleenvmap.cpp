#include <lightwave.hpp>

namespace lightwave {

class SimpleEnvMap final : public BackgroundLight {
    /// @brief The texture to use as background
    ref<Texture> m_texture;
    /// @brief An optional transform from local-to-world space
    ref<Transform> m_transform;

public:
    SimpleEnvMap(const Properties &properties)
    : BackgroundLight(properties) {
        m_texture   = properties.getChild<Texture>();
        m_transform = properties.getOptionalChild<Transform>();
    }

    DirectLightSample sampleDirect(const Point &origin,
                                   Sampler &rng) const override {
        PROFILE("Simple Direct")
        const Context cont = Context();
        auto sample_wi = squareToUniformSphere(rng.next2D());
        auto sample = sample_wi;
        if (m_transform) {
            SimpleTransform T = m_transform->interpolate(0);
            sample            = T.inverse(sample).normalized();
        }
        auto trans = Vector2((std::atan2(-sample.z(), sample.x()) + Pi) * Inv2Pi,
                    safe_acos(sample.y()) * InvPi);
        auto intensity = m_texture->evaluate(trans, cont);
        float pdf = 1.0 / (4.0 * Pi);
        return {
            .wi     = sample_wi,
            .weight = intensity / pdf,
            .pdf = pdf,
            .distance = Infinity,
        };
    }

    EmissionEval evaluate(const Vector &direction) const override {
        PROFILE("Simple Eval")
        const Context cont = Context();
        auto d = direction;
        if (m_transform) {
            SimpleTransform T = m_transform->interpolate(0);
            d                 = T.inverse(direction).normalized();
        }
        auto trans = Vector2((std::atan2(-d.z(), d.x()) + Pi) * Inv2Pi,
                            safe_acos(d.y()) * InvPi);
        auto intensity = m_texture->evaluate(trans, cont);
        float pdf = 1.0 / (4.0 * Pi);
        return {
            .value = intensity,
            .pdf = pdf,
        };
    }

    std::string toString() const override {
        return tfm::format(
            "SimpleEnvMap[\n"
            "   texture = %s,\n"
            "   transform = %s\n"
            "]",
            indent(m_texture),
            indent(m_transform));
    }
};

} // namespace lightwave

REGISTER_LIGHT(SimpleEnvMap, "simpleenvmap")
