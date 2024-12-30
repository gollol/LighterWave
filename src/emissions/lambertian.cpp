#include <lightwave.hpp>

namespace lightwave {

class Lambertian : public Emission {
    ref<Texture> m_emission;

public:
    Lambertian(const Properties &properties) {
        m_emission = properties.get<Texture>("emission");
    }

    EmissionEval evaluate(const Point2 &uv, const Vector &wo, Context &cont) const override {
        const float cosTerm = std::max(Frame::cosTheta(wo), 0.f);
        return {
            .value = cosTerm > 0 ? m_emission->evaluate(uv, cont) : Color::black(),
        };
    }

    std::string toString() const override {
        return tfm::format("Lambertian[\n"
                           "  emission = %s\n"
                           "]",
                           indent(m_emission));
    }
};

} // namespace lightwave

REGISTER_EMISSION(Lambertian, "lambertian")
