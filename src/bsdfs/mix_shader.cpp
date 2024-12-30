#include <lightwave.hpp>

namespace lightwave {

class MixShader : public Bsdf {
    ref<Texture> m_factor;
    ref<Bsdf> m_bsdf1;
    ref<Bsdf> m_bsdf2;

public:
    MixShader(const Properties &properties) {
        m_factor = properties.get<Texture>("factor");

        const std::vector<ref<Bsdf>> bsdf_children =
            properties.getChildren<Bsdf>();

        assert_condition(bsdf_children.size() == 2, {
            logger(EError,
                   "Exactly 2 BSDFs required for 'mix_shader' with ids: %s, %s, %s, %d given!",
                   this->id(),
                   m_factor->id(),
                   bsdf_children[0]->id(),
                   bsdf_children.size());
        });

        this->m_bsdf1 = bsdf_children[0];
        this->m_bsdf2 = bsdf_children[1];
    }

    Color albedo(const Point2 &uv, const Context &cont) const override {
        float factor = this->m_factor->scalar(uv, cont);

        assert_condition((factor >= 0.0f) and (factor <= 1.0f), {
            logger(EError,
                   "Factor must be a floating point number between 0 and 1, %f "
                   "was given!",
                   factor);
        });

        // Linearly interpolate value
        return factor * this->m_bsdf2->albedo(uv, cont) +
               (1.0f - factor) * this->m_bsdf2->albedo(uv, cont);
    }

    BsdfEval evaluate(const Point2 &uv, const Vector &wo, const Vector &wi,
                      const Context &cont, float roughness = 0) const override {
        BsdfEval bsdfEval1 = this->m_bsdf1->evaluate(uv, wo, wi, cont);
        BsdfEval bsdfEval2 = this->m_bsdf2->evaluate(uv, wo, wi, cont);
        float factor       = this->m_factor->scalar(uv, cont);

        assert_condition((factor >= 0.0f) and (factor <= 1.0f), {
            logger(EError,
                   "Factor must be a floating point number between 0 and 1, %f "
                   "was given!",
                   factor);
        });

        // Linearly interpolate both value and pdf of child BSDFs
        return BsdfEval{
            .value =
                factor * bsdfEval2.value + (1.0f - factor) * bsdfEval1.value,
            .pdf = factor * bsdfEval2.pdf + (1.0f - factor) * bsdfEval1.pdf,
        };
    }

    BsdfSample sample(const Point2 &uv, const Vector &wo, Sampler &rng,
                      const Context &cont, float roughness = 0) const override {
        float factor = this->m_factor->scalar(uv, cont);

        assert_condition((factor >= 0.0f) and (factor <= 1.0f), {
            logger(EError,
                   "Factor must be a floating point number between 0 and 1, %f "
                   "was given!",
                   factor);
        });

        if (rng.next() < factor) {
            // use second BSDF
            return this->m_bsdf2->sample(uv, wo, rng, cont);
        } else {
            // use first BSDF
            return this->m_bsdf1->sample(uv, wo, rng, cont);
        }
    }

    float getRoughness(const Point2 &uv, const Context &cont) const override {
        // TODO: check if this is correct
        return 0;
    }

    std::string toString() const override {
        return tfm::format(
            "MixShader[\n"
            "  factor   = %s,\n"
            "  BSDF1    = %s,\n"
            "  BSDF2    = %s\n"
            "]",
            indent(m_factor),
            indent(m_bsdf1),
            indent(m_bsdf2));
    }
};

} // namespace lightwave

REGISTER_BSDF(MixShader, "mix_shader")
