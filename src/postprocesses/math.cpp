#include <lightwave.hpp>

namespace lightwave {

class Math : public Postprocess {

    enum Variable {
        less,
        greater,
        multiply,
        divide,
        power,
    } m_operation;

    bool clamping;

    ref<Texture> b_image;

public:
    Math(const Properties &properties) : Postprocess(properties) {
        m_operation = properties.getEnum<Variable>("operation",
                                                   {
                                                       { "less", less },
                                                       { "greater", greater },
                                                       { "multiply", multiply },
                                                       { "divide", divide },
                                                       { "power", power },
                                                   });

        clamping = properties.get<bool>("clamping");
        b_image  = properties.get<Texture>("b_image");
    }

    /// @brief It is expected to always get a grey image where R=G=B and a R=G=B
    /// b_image. Does math...
    void execute() override {
        // TODO: check if this is correct
        Context cont;
        // Output initialized to resolution of input
        m_output->initialize(m_input->resolution());

        for (int y = 0; y < m_output->resolution().y(); y++)
            for (int x = 0; x < m_output->resolution().x(); x++) {
                Point2i pixel(x, y);
                Color rgb = m_input->get(pixel);
                float b_value =
                    b_image
                        ->evaluate(
                            Point2((float) x / m_output->resolution().x(),
                                   (float) (m_output->resolution().y() - y) /
                                       m_output->resolution().y()),
                            cont)
                        .r();

                Color out = Color(0);

                switch (m_operation) {
                case less:
                    out = Color(rgb.r() < b_value);
                    break;
                case greater:
                    out = Color(rgb.r() > b_value);
                    break;
                case multiply:
                    out = Color(rgb.r() * b_value);
                    break;
                case divide:
                    out = Color(rgb.r() / b_value);
                    break;
                case power:
                    out = Color(pow(rgb.r(), b_value));
                    break;
                }

                m_output->get(pixel) = out;
            }

        m_output->save();
    }

    std::string toString() const override {
        return tfm::format(
            "MathPP[\n"
            "  input = %s,\n"
            "  output = %s,\n"
            "  operation = %s, \n"
            "  clamping = %s, \n"
            "]",
            indent(m_input),
            indent(m_output),
            indent(m_operation),
            indent(clamping)
        );
    }
};

} // namespace lightwave

REGISTER_POSTPROCESS(Math, "math");
