#include <lightwave.hpp>

namespace lightwave
{

    class SunLight final : public BackgroundLight
    {
        Vector m_direction;
        Color m_intensity;
        float m_size;
        float pdf;

    public:
        SunLight(const Properties &properties)
            : BackgroundLight(properties)
        {
            m_direction = properties.get<Vector>("direction").normalized();
            m_intensity = properties.get<Color>("intensity");
            m_size = properties.get<float>("size");
            pdf = uniformPartialHemispherePdfConst(m_size);
        }

        DirectLightSample sampleDirect(const Point &origin,
                                       Sampler &rng) const override
        {
            auto sample = squareToUniformPartialHemisphere(rng.next2D(), m_size);
            auto transformed_sample = transformFromZ(sample, m_direction);
            return {
                .wi = transformed_sample,
                .weight = m_intensity,
                .pdf = pdf,
                .distance = Infinity,
            };
        }

        EmissionEval evaluate(const Vector &direction) const override
        {
            float cos_angle = direction.dot(m_direction);
            float angle = acos(cos_angle);
            float current_pdf = 0;
            Color current_value(0.0, 0.0, 0.0);
            if (angle <= m_size)
            {
                current_pdf = pdf;
                current_value = m_intensity * pdf;
            }
            return {
                .value = current_value,
                .pdf = current_pdf,
            };
        }

        std::string toString() const override
        {
            return tfm::format("SunLight[\n"
                               "]");
        }
    };

} // namespace lightwave

REGISTER_LIGHT(SunLight, "sun")
