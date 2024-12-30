#include <lightwave.hpp>

namespace lightwave {

class Fog : public Postprocess {
    //TODO set fog color
    ref<Image> distance;

public:
    Fog(const Properties &properties) : Postprocess(properties) {
        NOT_IMPLEMENTED
        distance = properties.get<Image>("distance");
    }

    /// @brief Tonemapping of image...
    void execute() override {
        // Output initialized to resolution of input
        m_output->initialize(m_input->resolution());

        // get max rendered distance
        float d_max = -Infinity;
        for (int y = 0; y < m_output->resolution().y(); y++)
            for (int x = 0; x < m_output->resolution().x(); x++) {
                Point2i pixel(x, y);
                float d = distance->get(pixel).r();

                if (std::isfinite(d))
                    if (d_max < d)
                        d_max = d;
            }

        // add fog to scene with linerar interpolation of a fog "image" and the
        // m_input
        for (int y = 0; y < m_output->resolution().y(); y++)
            for (int x = 0; x < m_output->resolution().x(); x++) {
                Point2i pixel(x, y);
                Color rgb = m_input->get(pixel);

                float d = distance->get(pixel).r();

                // Only apply fog where rays hit the scene
                if (std::isfinite(d)) {
                    float fog = 3 * powf(d / d_max, 1);
                    fog       = clamp(fog, 0.0f, 1.0f);

                    rgb = (1.0f - fog) * rgb + fog * Color(1);
                }

                m_output->get(pixel) = rgb;
            }

        m_output->save();
    }

    std::string toString() const override {
        return tfm::format(
            "FogPP[\n"
            "  input = %s,\n"
            "  output = %s,\n"
            "]",
            indent(m_input),
            indent(m_output)
        );
    }
};

} // namespace lightwave

REGISTER_POSTPROCESS(Fog, "fog");
