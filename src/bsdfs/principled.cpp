#include <lightwave.hpp>

#include "fresnel.hpp"
#include "microfacet.hpp"

namespace lightwave {

/**
 * Lobes required:
 *  - Diffuse
 *  - Conductor
 *  - Transmission
 *  - (Specular)
 */

/// @brief Calculates the reflection coefficient for light incoming parallel to
/// the normal (θ = 0)
/// @param eta The index of refraction for the second material, assuming the
/// first is air
/// @return The reflection coefficient for light incoming parallel to the normal
float F0_from_ior(const float eta) {
    float f0 = (eta - 1.0f) / (eta + 1.0f);
    return sqr(f0);
}

float ior_from_F0(float f0) {
    float sqrt_f0 = sqrt(clamp(f0, 0.f, 0.99f));
    return (1.0f + sqrt_f0) / (1.0f - sqrt_f0);
}

struct DiffuseLobe {
    Color color;

    BsdfEval evaluate(const Vector &wo, const Vector &wi) const {
        if (!Frame::sameHemisphere(wo, wi))
            return BsdfEval::invalid();
        const float cosTerm = Frame::absCosTheta(wi);
        return {
            .value = color * (cosTerm * InvPi),
            .pdf   = cosTerm * InvPi,
        };
        // hints:
        // * copy your diffuse bsdf evaluate here
        // * you do not need to query a texture, the albedo is given by `color`
    }

    BsdfSample sample(const Vector &wo, Sampler &rng) const {
        const Vector wi = squareToCosineHemisphere(rng.next2D());

        return {
            .wi     = wi * (Frame::cosTheta(wi) > 0.0f ? +1.0f : -1.0f),
            .weight = color,
            .pdf    = Frame::cosTheta(wi) * InvPi,
        };
        // hints:
        // * copy your diffuse bsdf evaluate here
        // * you do not need to query a texture, the albedo is given by `color`
    }
};

struct SpecularLobe {
    float alpha;
    Color color;
    float f0;

    BsdfEval evaluate(const Vector &wo, const Vector &wi) const {
        const auto normal = (wi + wo).normalized();

        // VNDF PDF
        float pdf = microfacet::pdfGGXVNDF(alpha, normal, wo);
        if (!(pdf > 0))
            return BsdfEval::invalid();
        pdf *= microfacet::detReflection(normal, wo);

        const float Gi =
            microfacet::anisotropicSmithG1(alpha, alpha, normal, wi);
        return {
            .value = color * (Gi * pdf) * schlick(f0, abs(normal.dot(wo))),
            .pdf   = pdf,
        };
        // hints:
        // * copy your roughconductor bsdf evaluate here
        // * you do not need to query textures
        //   * the reflectance is given by `color'
        //   * the variable `alpha' is already provided for you
    }

    BsdfSample sample(const Vector &wo, Sampler &rng) const {
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
            .weight = color * Gi * schlick(f0, abs(normal.dot(wo))),
            .pdf    = pdf,
        };
        // hints:
        // * copy your roughconductor bsdf sample here
        // * you do not need to query textures
        //   * the reflectance is given by `color'
        //   * the variable `alpha' is already provided for you
    }
};

struct TransmissionLobe {
    float alpha;
    Color color;
    float ior;

    BsdfEval evaluate(const Vector &wo, const Vector &wi) const {
        if (Frame::sameHemisphere(wi, wo)) {
            // We're not in the same hemisphere, so we can't evaluate
            // transmission
            return BsdfEval::invalid();
        }

        const float cosTheta = Frame::cosTheta(wo);
        float eta            = cosTheta < 0 ? 1 / ior : ior;
        const Vector normal  = (wi * eta + wo).normalized();

        // VNDF PDF
        float pdf = microfacet::pdfGGXVNDF(alpha, normal, wo);

        if (!(pdf > 0)) {
            return BsdfEval::invalid();
        }

        float Gi = microfacet::smithG1(alpha, normal, wi);
        pdf *= microfacet::detRefraction(normal, wi, wo, eta);

        return {
            .value = color * Color(pdf * Gi / sqr(eta)),
            .pdf   = pdf,
        };
    }

    BsdfSample sample(const Vector &wo, Sampler &rng) const {
        const float cosTheta = Frame::cosTheta(wo);
        float eta            = cosTheta < 0 ? 1 / ior : ior;

        const Vector normal =
            microfacet::sampleGGXVNDF(alpha, wo, rng.next2D());
        float pdf       = microfacet::pdfGGXVNDF(alpha, normal, wo);
        const Vector wi = refract(wo, normal, eta);

        if (Frame::sameHemisphere(wi, wo)) {
            return BsdfSample::invalid();
        }

        const float Gi =
            microfacet::anisotropicSmithG1(alpha, alpha, normal, wi);

        return {
            .wi     = wi,
            .weight = color * Color(Gi / sqr(eta)),
            .pdf    = pdf * microfacet::detRefraction(normal, wi, wo, eta),
        };
    }
};

class Principled : public Bsdf {
    ref<Texture> m_baseColor;
    ref<Texture> m_roughness;
    ref<Texture> m_metallic;
    ref<Texture> m_ior;
    ref<Texture> m_specularLevel;
    ref<Texture> m_specularTint;
    ref<Texture> m_transmission;

    struct Combination {
        float diffuseSelectionProb;
        float metallicSelectionProb;
        float transmissionSelectionProb;
        DiffuseLobe diffuse;
        SpecularLobe specular;
        TransmissionLobe transmission;
    };

    Combination combine(const Point2 &uv, const Vector &wo, const Context &cont,
                        float roughness) const {
        const auto baseColor = m_baseColor->evaluate(uv, cont);
        const auto alpha     = max(
            float(1e-3), sqr(max(m_roughness->scalar(uv, cont), roughness)));
        const auto specularLevel = m_specularLevel->scalar(uv, cont);
        const auto metallic      = m_metallic->scalar(uv, cont);
        const auto transmission  = m_transmission->scalar(uv, cont);
        const auto specularTint  = m_specularTint->evaluate(uv, cont);
        const float ior = std::max(float(1e-5), m_ior->scalar(uv, cont));

        SpecularLobe specularLobe = {
            .alpha = alpha,
            .color = Color(0.0f),
            .f0    = 0.0f,
        };
        TransmissionLobe transmissionLobe = {
            .alpha = alpha,
            .color = Color(0.0f),
            .ior   = ior,
        };
        DiffuseLobe diffuseLobe = {
            .color = Color(0.0f),
        };

        Color weight = Color(1.0f); // Current weight to be used and decreased
                                    // for each lobe

        // Metallic component
        specularLobe.color = weight * metallic;

        weight *= std::max((1.0f - metallic), 0.0f);

        // Transmission component
        specularLobe.color +=
            weight * transmission * F0_from_ior(ior) * specularTint;
        transmissionLobe.color =
            transmission * (Color(1.0f) - (F0_from_ior(ior) * specularTint)) *
            baseColor;

        weight *= std::max((1.0f - transmission), 0.0f);

        // Specular component
        float eta = ior;
        float f0  = F0_from_ior(eta);
        if (specularLevel != 0.5f) {
            f0 *= 2.0f * specularLevel;
            eta = ior_from_F0(f0);
            if (ior < 1.0f) {
                eta = 1.0f / eta;
            }
        }

        specularLobe.f0 = f0;

        specularLobe.color += specularTint;
        // weight *= std::max((1.0f - f0), 0.0f);

        // Diffuse component
        diffuseLobe.color = weight * baseColor;

        const auto diffuseAlbedo      = diffuseLobe.color.mean();
        const auto metallicAlbedo     = specularLobe.color.mean();
        const auto transmissionAlbedo = transmissionLobe.color.mean();

        const auto totalAlbedo =
            diffuseAlbedo + metallicAlbedo + transmissionAlbedo;

        return {
            .diffuseSelectionProb =
                totalAlbedo > 0 ? diffuseAlbedo / totalAlbedo : 1.0f,
            .metallicSelectionProb =
                totalAlbedo > 0 ? metallicAlbedo / totalAlbedo : 0.0f,
            .transmissionSelectionProb =
                totalAlbedo > 0 ? transmissionAlbedo / totalAlbedo : 0.0f,
            .diffuse      = diffuseLobe,
            .specular     = specularLobe,
            .transmission = transmissionLobe,
        };
    }

public:
    Principled(const Properties &properties) {
        m_baseColor     = properties.get<Texture>("baseColor");
        m_metallic      = properties.get<Texture>("metallic");
        m_roughness     = properties.get<Texture>("roughness");
        m_ior           = properties.get<Texture>("ior");
        m_specularLevel = properties.get<Texture>("specularLevel");
        m_specularTint  = properties.get<Texture>("specularTint");
        m_transmission  = properties.get<Texture>("transmission");
    }

    Color albedo(const Point2 &uv, const Context &cont) const override {
        return m_baseColor->evaluate(uv, cont);
    }

    BsdfEval evaluate(const Point2 &uv, const Vector &wo, const Vector &wi,
                      const Context &cont, float roughness = 0) const override {
        PROFILE("Principled")

        const auto combination  = combine(uv, wo, cont, roughness);
        const auto diffuse      = combination.diffuse.evaluate(wo, wi);
        const auto specular     = combination.specular.evaluate(wo, wi);
        const auto transmission = combination.transmission.evaluate(wo, wi);

        const auto pdf =
            diffuse.pdf * combination.diffuseSelectionProb +
            specular.pdf * combination.metallicSelectionProb +
            transmission.pdf * combination.transmissionSelectionProb;

        if (!(pdf > 0))
            return BsdfEval::invalid();

        return {
            .value = diffuse.value + specular.value + transmission.value,
            .pdf   = pdf,
        };
        // hint: evaluate `combination.diffuse` and `combination.metallic` and
        // combine their results
    }

    BsdfSample sample(const Point2 &uv, const Vector &wo, Sampler &rng,
                      const Context &cont, float roughness = 0) const override {
        PROFILE("Principled")

        const auto combination = combine(uv, wo, cont, roughness);

        const auto rngVal = rng.next();

        BsdfSample sample;

        if (rngVal < combination.diffuseSelectionProb) {
            sample = combination.diffuse.sample(wo, rng);
        } else if ((rngVal >= combination.diffuseSelectionProb) and
                   (rngVal < (combination.diffuseSelectionProb +
                              combination.metallicSelectionProb))) {
            sample = combination.specular.sample(wo, rng);
        } else {
            sample = combination.transmission.sample(wo, rng);
        }

        if (sample.isInvalid())
            return sample;

        const auto diffuse  = combination.diffuse.evaluate(wo, sample.wi);
        const auto specular = combination.specular.evaluate(wo, sample.wi);
        const auto transmission =
            combination.transmission.evaluate(wo, sample.wi);

        const auto pdf =
            diffuse.pdf * combination.diffuseSelectionProb +
            specular.pdf * combination.metallicSelectionProb +
            transmission.pdf * combination.transmissionSelectionProb;

        if (!(pdf > 0)) {
            return BsdfSample::invalid();
        }


        return {
            .wi = sample.wi,
            .weight =
                (diffuse.value + specular.value + transmission.value) / pdf,
            .pdf = pdf,
        };
        // hint: sample either `combination.diffuse` (probability
        // `combination.diffuseSelectionProb`) or `combination.metallic`
    }

    float getRoughness(const Point2 &uv, const Context &cont) const override {
        return m_roughness->scalar(uv, cont);
    }

    std::string toString() const override {
        return tfm::format(
            "Principled[\n"
            "  baseColor     = %s,\n"
            "  metallic      = %s,\n"
            "  roughness     = %s,\n"
            "  ior           = %s,\n"
            "  specularLevel = %s,\n"
            "  specularTint  = %s,\n"
            "  transmission  = %s,\n"
            "]",
            indent(m_baseColor),
            indent(m_metallic),
            indent(m_roughness),
            indent(m_ior),
            indent(m_specularLevel),
            indent(m_specularTint),
            indent(m_transmission));
    }
};

} // namespace lightwave

REGISTER_BSDF(Principled, "principled")
