#include "fresnel.hpp"
#include <lightwave.hpp>

namespace lightwave {

class Dielectric : public Bsdf {
    ref<Texture> m_ior;
    ref<Texture> m_reflectance;
    ref<Texture> m_transmittance;

public:
    Dielectric(const Properties &properties) {
        m_ior           = properties.get<Texture>("ior");
        m_reflectance   = properties.get<Texture>("reflectance");
        m_transmittance = properties.get<Texture>("transmittance");
    }

    Color albedo(const Point2 &uv, const Context &cont) const override {
        return m_reflectance->evaluate(uv, cont) +
               m_transmittance->evaluate(uv, cont);
    }

    BsdfEval evaluate(const Point2 &uv, const Vector &wo, const Vector &wi,
                      const Context &cont, float roughness = 0) const override {
        // the probability of a light sample picking exactly the direction `wi'
        // that results from reflecting or refracting `wo' is zero, hence we can
        // just ignore that case and always return black
        return BsdfEval::invalid();
    }

    BsdfSample sample(const Point2 &uv, const Vector &wo, Sampler &rng,
                      const Context &cont, float roughness = 0) const override {
        BsdfSample result;

        float cosTheta = Frame::cosTheta(wo);
        float eta      = m_ior->scalar(uv, cont);

        if (cosTheta < 0) {
            eta = 1 / eta;
        }

        const float F        = fresnelDielectric(cosTheta, eta);
        const auto cReflect  = F * m_reflectance->evaluate(uv, cont);
        const auto cTransmit = (1 - F) * m_transmittance->evaluate(uv, cont);
        const float pReflect =
            cReflect.mean() > 0
                ? cReflect.mean() / (cReflect.mean() + cTransmit.mean())
                : 0;
        if (rng.next() < pReflect) {
            // sample reflection
            result.wi     = reflect(wo, Vector(0, 0, 1));
            result.weight = cReflect / pReflect;
        } else {
            // sample refraction
            result.wi     = refract(wo, Vector(0, 0, 1), eta);
            result.weight = cTransmit / ((1 - pReflect) * (eta * eta));

            if (result.wi.isZero())
                return BsdfSample::invalid();
        }
        result.pdf = Infinity;
        assert_finite(result.wi, {});
        return result;
    }

    float getRoughness(const Point2 &uv, const Context &cont) const override {
        return 0;
    }

    std::string toString() const override {
        return tfm::format(
            "Dielectric[\n"
            "  ior           = %s,\n"
            "  reflectance   = %s,\n"
            "  transmittance = %s\n"
            "]",
            indent(m_ior),
            indent(m_reflectance),
            indent(m_transmittance));
    }
};

} // namespace lightwave

REGISTER_BSDF(Dielectric, "dielectric")
