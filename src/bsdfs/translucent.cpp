#include <lightwave.hpp>

namespace lightwave {

class Translucent : public Bsdf {
    ref<Texture> m_albedo;

public:
    Translucent(const Properties &properties) {
        m_albedo = properties.get<Texture>("albedo");
    }

    Color albedo(const Point2 &uv, const Context &cont) const override {
        return m_albedo->evaluate(uv, cont);
    }

    BsdfEval evaluate(const Point2 &uv, const Vector &wo, const Vector &wi,
                      const Context &cont, float roughness = 0) const override {
        Vector wii = wi;
        if (!Frame::sameHemisphere(wo, wi))
            wii = -wi;
        const float cosTerm = Frame::absCosTheta(wii);
        return {
            .value = cosTerm * m_albedo->evaluate(uv, cont) * InvPi,
            .pdf   = cosTerm * InvPi,
        };
    }

    BsdfSample sample(const Point2 &uv, const Vector &wo, Sampler &rng,
                      const Context &cont, float roughness = 0) const override {
        Vector wi = squareToCosineHemisphere(rng.next2D());
        return {
            .wi     = wi * (Frame::cosTheta(wo) > 0.0f ? -1.0f : +1.0f),
            .weight = m_albedo->evaluate(uv, cont),
            .pdf    = Frame::cosTheta(wi) * InvPi,
        };
    }

    float getRoughness(const Point2 &uv, const Context &cont) const override {
        return 1;
    }

    std::string toString() const override {
        return tfm::format(
            "Translucent[\n"
            "  albedo = %s\n"
            "]",
            indent(m_albedo));
    }
};

} // namespace lightwave

REGISTER_BSDF(Translucent, "translucent")
