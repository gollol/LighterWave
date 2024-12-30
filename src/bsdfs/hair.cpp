#include <lightwave.hpp>

namespace lightwave {

enum Component { REFLECTION, TRANSMISSION };

/// @brief Currently just a simple Diffuse
class Hair : public Bsdf {

    ref<Texture> m_color;
    ref<Texture> m_offset;
    ref<Texture> m_roughnessU;
    ref<Texture> m_roughnessV;
    ref<Texture> m_tangent;
    Component m_component;

public:
    Hair(const Properties &properties) {

        m_color      = properties.get<Texture>("color");
        m_offset     = properties.get<Texture>("offset");
        m_roughnessU = properties.get<Texture>("roughnessU");
        m_roughnessV = properties.get<Texture>("roughnessV");
        m_tangent    = properties.get<Texture>("tangent");
        m_component  = properties.getEnum<Component>(
            "component",
            { { "reflection", REFLECTION }, { "transmission", TRANSMISSION } });
    }

    Color albedo(const Point2 &uv, const Context &cont) const override {
        return m_color->evaluate(uv, cont);
    }

    BsdfEval evaluate(const Point2 &uv, const Vector &wo, const Vector &wi,
                      const Context &cont, float roughness = 0) const override {
        PROFILE("Hair")
        if (!Frame::sameHemisphere(wo, wi))
            return BsdfEval::invalid();
        const float cosTerm = Frame::absCosTheta(wi);
        return {
            .value = cosTerm * m_color->evaluate(uv, cont) * InvPi,
            .pdf   = cosTerm * InvPi,
        };
    }

    BsdfSample sample(const Point2 &uv, const Vector &wo, Sampler &rng,
                      const Context &cont, float roughness = 0) const override {
        const Vector wi = squareToCosineHemisphere(rng.next2D());
        return {
            .wi     = wi * (Frame::cosTheta(wo) > 0.0f ? +1.0f : -1.0f),
            .weight = m_color->evaluate(uv, cont),
            .pdf    = Frame::cosTheta(wi) * InvPi,
        };
    }

    float getRoughness(const Point2 &uv, const Context &cont) const override {
        // TODO: check if this is correct
        return 0;
    }

    std::string toString() const override {
        return tfm::format(
            "Hair[\n"
            "  color = %s\n"
            "  offset = %s\n"
            "  roughnessU = %s\n"
            "  roughnessV = %s\n"
            "  tangent = %s\n"
            "  component = %s\n"
            "]",
            indent(m_color),
            indent(m_offset),
            indent(m_roughnessU),
            indent(m_roughnessV),
            indent(m_tangent),
            indent(m_component));
    }
};

} // namespace lightwave

REGISTER_BSDF(Hair, "hair")
