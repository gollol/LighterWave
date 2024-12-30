#include <lightwave.hpp>

namespace lightwave {

class AddShader : public Bsdf {
    ref<Bsdf> m_bsdf1;
    ref<Bsdf> m_bsdf2;

public:
    AddShader(const Properties &properties) {
        const std::vector<ref<Bsdf>> bsdf_children =
            properties.getChildren<Bsdf>();

        assert_condition(bsdf_children.size() == 2, {
            logger(EError,
                   "Exactly 2 BSDFs required for 'add_shader', %d given!",
                   bsdf_children.size());
        });

        this->m_bsdf1 = bsdf_children[0];
        this->m_bsdf2 = bsdf_children[1];
    }

    Color albedo(const Point2 &uv, const Context &cont) const override {
        return this->m_bsdf2->albedo(uv, cont) +
               this->m_bsdf2->albedo(uv, cont);
    }

    BsdfEval evaluate(const Point2 &uv, const Vector &wo, const Vector &wi,
                      const Context &cont, float roughness = 0) const override {
        BsdfEval bsdfEval1 = this->m_bsdf1->evaluate(uv, wo, wi, cont);
        BsdfEval bsdfEval2 = this->m_bsdf2->evaluate(uv, wo, wi, cont);

        return BsdfEval{
            .value = bsdfEval2.value + bsdfEval1.value,
            .pdf   = (bsdfEval2.pdf + bsdfEval1.pdf) * 0.5f,
        };
    }

    BsdfSample sample(const Point2 &uv, const Vector &wo, Sampler &rng,
                      const Context &cont, float roughness = 0) const override {
        auto out_shad = rng.next() < 0.5
                            ? this->m_bsdf2->sample(uv, wo, rng, cont)
                            : this->m_bsdf1->sample(uv, wo, rng, cont);
        out_shad.pdf *= 0.5f;

        // auto E = evaluate(uv, wo, out_shad.wi);
        // return {
        //     .weight = E.value / E.pdf,
        //     .pdf = E.pdf
        // };

        return out_shad;
    }

    float getRoughness(const Point2 &uv, const Context &cont) const override {
        // TODO: check if this is correct
        return 0;
    }

    std::string toString() const override {
        return tfm::format(
            "AddShader[\n"
            "  BSDF1    = %s,\n"
            "  BSDF2    = %s\n"
            "]",
            indent(m_bsdf1),
            indent(m_bsdf2));
    }
};

} // namespace lightwave

REGISTER_BSDF(AddShader, "add_shader")
