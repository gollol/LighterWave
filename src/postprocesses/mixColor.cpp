#include <lightwave.hpp>

namespace lightwave {

class MixColor : public Postprocess {
    ref<Texture> b_image;

    enum Variable {
        mix,
        add,
        multiply,
    } m_operation;
    float factor;
    bool clamping;

public:
    MixColor(const Properties &properties) : Postprocess(properties) {
        m_operation = properties.getEnum<Variable>("operation",
                                                   {
                                                       { "mix", mix },
                                                       { "add", add },
                                                       { "multiply", multiply },
                                                   });

        factor   = properties.get<float>("factor");
        clamping = properties.get<bool>("clamping");

        b_image = properties.get<Texture>("b_image");
    }

    /// @brief It combines images/colors together dependent of the operation.
    /// m_image should be OG image. B_Image gets mixed with OG image.
    void execute() override {
        // TODO: check if this is correct
        Context cont;
        // Output initialized to resolution of input
        m_output->initialize(m_input->resolution());

        for (int y = 0; y < m_output->resolution().y(); y++)
            for (int x = 0; x < m_output->resolution().x(); x++) {
                Point2i pixel(x, y);
                Color rgb = m_input->get(pixel);

                Color out   = Color(0);
                Color b_rgb = b_image->evaluate(
                    Point2((float) x / m_output->resolution().x(),
                           (float) (m_output->resolution().y() - y) /
                               m_output->resolution().y()),
                    cont);

                switch (m_operation) {
                case mix:
                    out = (1.0f - factor) * rgb + factor * b_rgb;
                    break;
                case add:
                    out = rgb + factor * b_rgb;
                    break;
                case multiply:
                    Color multiplied = rgb * b_rgb;
                    out = (1.0f - factor) * rgb + factor * multiplied;
                    break;
                }

                if (clamping) {
                    out.r() = clamp(out.r(), 0.0f, 1.0f);
                    out.g() = clamp(out.g(), 0.0f, 1.0f);
                    out.b() = clamp(out.b(), 0.0f, 1.0f);
                }

                m_output->get(pixel) = out;
            }

        m_output->save();
    }

    std::string toString() const override { return "Color mixed, BAMM!"; }
};

} // namespace lightwave

REGISTER_POSTPROCESS(MixColor, "mixColor");
