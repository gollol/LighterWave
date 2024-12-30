#include "fresnel.hpp"
#include "microfacet.hpp"
#include <lightwave.hpp>

namespace lightwave {

class RoughConductor : public Bsdf {
    ref<Texture> m_reflectance;
    ref<Texture> m_roughness;

public:
    RoughConductor(const Properties &properties) {
        m_reflectance = properties.get<Texture>("reflectance");
        m_roughness   = properties.get<Texture>("roughness");
    }

    Color albedo(const Point2 &uv, const Context &cont) const override {
        return m_reflectance->evaluate(uv, cont);
    }

    BsdfEval evaluate(const Point2 &uv, const Vector &wo, const Vector &wi,
                      const Context &cont, float roughness = 0) const override {
        // Using the squared roughness parameter results in a more gradual
        // transition from specular to rough. For numerical stability, we avoid
        // extremely specular distributions (alpha values below 10^-3)
        const auto alpha = std::max(
            float(1e-3), sqr(max(m_roughness->scalar(uv, cont), roughness)));

        const auto normal = (wi + wo).normalized();

        // VNDF PDF
        float pdf = microfacet::pdfGGXVNDF(alpha, normal, wo);
        if (!(pdf > 0))
            return BsdfEval::invalid();
        pdf *= microfacet::detReflection(normal, wo);

        const float Gi =
            microfacet::anisotropicSmithG1(alpha, alpha, normal, wi);
        return {
            .value = m_reflectance->evaluate(uv, cont) * (Gi * pdf),
            .pdf   = pdf,
        };
        // hints:
        // * the microfacet normal can be computed from `wi' and `wo'
    }

    BsdfSample sample(const Point2 &uv, const Vector &wo, Sampler &rng,
                      const Context &cont, float roughness = 0) const override {
        const auto alpha = std::max(
            float(1e-3), sqr(max(m_roughness->scalar(uv, cont), roughness)));

        const Vector normal =
            microfacet::sampleGGXVNDF(alpha, wo, rng.next2D());

        const Vector wi = reflect(wo, normal);
        if (!Frame::sameHemisphere(wi, wo))
            return BsdfSample::invalid();

        // VNDF PDF
        float pdf = microfacet::pdfGGXVNDF(alpha, normal, wo);
        if (!(pdf > 0))
            return BsdfSample::invalid();
        pdf *= microfacet::detReflection(normal, wo);

        const float Gi = microfacet::smithG1(alpha, normal, wi);
        return {
            .wi     = wi,
            .weight = m_reflectance->evaluate(uv, cont) * Gi,
            .pdf    = pdf,
        };
        // hints:
        // * do not forget to cancel out as many terms from your equations as
        // possible!
        //   (the resulting sample weight is only a product of two factors)
    }

    float getRoughness(const Point2 &uv, const Context &cont) const override {
        return m_roughness->scalar(uv, cont);
    }

    std::string toString() const override {
        return tfm::format(
            "RoughConductor[\n"
            "  reflectance = %s,\n"
            "  roughness = %s\n"
            "]",
            indent(m_reflectance),
            indent(m_roughness));
    }
};

} // namespace lightwave

REGISTER_BSDF(RoughConductor, "roughconductor")
