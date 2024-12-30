#include <lightwave.hpp>

namespace lightwave {

class Conductor : public Bsdf {
    ref<Texture> m_reflectance;

public:
    Conductor(const Properties &properties) {
        m_reflectance = properties.get<Texture>("reflectance");
    }

    Color albedo(const Point2 &uv, const Context &cont) const override {
        return m_reflectance->evaluate(uv, cont);
    }

    BsdfEval evaluate(const Point2 &uv, const Vector &wo, const Vector &wi,
                      const Context &cont, float roughness = 0) const override {
        // the probability of a light sample picking exactly the direction `wi'
        // that results from reflecting `wo' is zero, hence we can just ignore
        // that case and always return black
        return BsdfEval::invalid();
    }

    BsdfSample sample(const Point2 &uv, const Vector &wo, Sampler &rng,
                      const Context &cont, float roughness = 0) const override {
        return {
            .wi     = reflect(wo, Vector(0, 0, 1)),
            .weight = m_reflectance->evaluate(uv, cont),
            .pdf    = Infinity,
        };
    }

    float getRoughness(const Point2 &uv, const Context &cont) const override {
        return 0;
    }

    std::string toString() const override {
        return tfm::format(
            "Conductor[\n"
            "  reflectance = %s\n"
            "]",
            indent(m_reflectance));
    }
};

} // namespace lightwave

REGISTER_BSDF(Conductor, "conductor")
