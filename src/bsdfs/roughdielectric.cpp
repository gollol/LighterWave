#include <lightwave.hpp>

#include "fresnel.hpp"
#include "microfacet.hpp"

namespace lightwave {

class RoughDielectric : public Bsdf {
    ref<Texture> m_ior;
    ref<Texture> m_roughness;
    ref<Texture> m_reflectance;
    ref<Texture> m_transmittance;

public:
    RoughDielectric(const Properties &properties) {
        m_ior           = properties.get<Texture>("ior");
        m_roughness     = properties.get<Texture>("roughness");
        m_reflectance   = properties.get<Texture>("reflectance");
        m_transmittance = properties.get<Texture>("transmittance");
    }

    Color albedo(const Point2 &uv, const Context &cont) const override {
        return m_reflectance->evaluate(uv, cont) +
               m_transmittance->evaluate(uv, cont);
    }

    BsdfEval evaluate(const Point2 &uv, const Vector &wo, const Vector &wi,
                      const Context &cont, float roughness = 0) const override {
        const bool isReflection = Frame::sameHemisphere(wi, wo);
        const float cosTheta    = Frame::cosTheta(wo);
        float eta               = m_ior->scalar(uv, cont);
        if (cosTheta < 0) {
            eta = 1 / eta;
        }

        const auto alpha = std::max(
            float(1e-3), sqr(max(m_roughness->scalar(uv, cont), roughness)));
        const Vector normal = isReflection ? (wi + wo).normalized()
                                           : (wi * eta + wo).normalized();

        // VNDF PDF
        float pdf = microfacet::pdfGGXVNDF(alpha, normal, wo);
        if (!(pdf > 0))
            return BsdfEval::invalid();

        const float Gi = microfacet::smithG1(alpha, normal, wi);
        const float F  = fresnelDielectric(normal.dot(wo), eta);
        if (isReflection) {
            pdf *= F * microfacet::detReflection(normal, wo);
            return {
                .value = Color(pdf * Gi) * m_reflectance->evaluate(uv, cont),
                .pdf   = pdf,
            };
        } else {
            pdf *= (1 - F) * microfacet::detRefraction(normal, wi, wo, eta);
            return {
                .value = Color(pdf * Gi / sqr(eta)) *
                         m_transmittance->evaluate(uv, cont),
                .pdf = pdf,
            };
        }
    }

    BsdfSample sample(const Point2 &uv, const Vector &wo, Sampler &rng,
                      const Context &cont, float roughness = 0) const override {
        const float cosTheta = Frame::cosTheta(wo);
        float eta            = m_ior->scalar(uv, cont);
        if (cosTheta < 0) {
            eta = 1 / eta;
        }

        const auto alpha = std::max(
            float(1e-3), sqr(max(m_roughness->scalar(uv, cont), roughness)));
        const Vector normal =
            microfacet::sampleGGXVNDF(alpha, wo, rng.next2D());

        float pdf = microfacet::pdfGGXVNDF(alpha, normal, wo);
        if (!(pdf > 0))
            return BsdfSample::invalid();

        const float F        = fresnelDielectric(normal.dot(wo), eta);
        const auto cReflect  = F * m_reflectance->evaluate(uv, cont);
        const auto cTransmit = (1 - F) * m_transmittance->evaluate(uv, cont);
        const float pReflect =
            cReflect.mean() > 0
                ? cReflect.mean() / (cReflect.mean() + cTransmit.mean())
                : 0;
        if (rng.next() < pReflect) {
            // sample reflection
            const Vector wi = reflect(wo, normal).normalized();
            if (!Frame::sameHemisphere(wi, wo))
                return BsdfSample::invalid();
            const float Gi =
                microfacet::anisotropicSmithG1(alpha, alpha, normal, wi);
            pdf *= F * microfacet::detReflection(normal, wo); 
            return {
                .wi     = wi,
                .weight = Color(Gi) * cReflect / pReflect,
                .pdf    = pdf,
            };
        } else {
            // sample transmission
            const Vector wi = refract(wo, normal, eta);
            if (wi.isZero() || Frame::sameHemisphere(wi, wo))
                return BsdfSample::invalid();
                
            const float Gi =
                microfacet::anisotropicSmithG1(alpha, alpha, normal, wi);
                 
            pdf *= (1 - F) * microfacet::detRefraction(normal, wi, wo, eta);
             
            return {
                .wi     = wi,
                .weight = Color(Gi / sqr(eta)) * cTransmit / (1 - pReflect),
                .pdf    = pdf,
            };
        }
    }

    float getRoughness(const Point2 &uv, const Context &cont) const override {
        return m_roughness->scalar(uv, cont);
    }

    std::string toString() const override {
        return tfm::format(
            "RoughDielectric[\n"
            "  ior           = %s,\n"
            "  roughness     = %s,\n"
            "  reflectance   = %s,\n"
            "  transmittance = %s\n"
            "]",
            indent(m_ior),
            indent(m_roughness),
            indent(m_reflectance),
            indent(m_transmittance));
    }
};

} // namespace lightwave

REGISTER_BSDF(RoughDielectric, "roughdielectric")
