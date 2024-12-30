#include <lightwave.hpp>

namespace lightwave {

class Color_Balance : public Postprocess {
    ref<Texture> whitepoint;

public:
    Color_Balance(const Properties &properties) : Postprocess(properties) {
        whitepoint = properties.get<Texture>("whitepoint");
    }

    /// @brief Image gets more/less contrast/brightness
    void execute() override {
        // TODO: check if this is correct
        Context cont;
        // Output initialized to resolution of input
        m_output->initialize(m_input->resolution());

        // A way to somehow calculate a Whitepoint
        // Color sum = Color(0);
        //
        // Loop over each Pixel
        // for(int y = 0; y < m_output->resolution().y(); y++)
        //{
        //     for(int x = 0; x < m_output->resolution().x(); x++)
        //     {
        //         Point2i pixel(x,y);
        //         Color rgb = m_input->get(pixel);
        //
        //         float maxIntensity = 0.0f;
        //         bool isEmissive = false;
        //         if(rgb.r() > 1.0f || rgb.g() > 1.0f || rgb.b() > 1.0f)
        //         {
        //             maxIntensity = max(max(rgb.r(), rgb.g()), rgb.b());
        //             isEmissive = true;
        //         }
        //
        //         if(isEmissive)
        //             rgb /= (maxIntensity + 0.001f);
        //
        //         sum += rgb;
        //
        //         //m_output->get(pixel) = rgb;
        //     }
        // }
        //
        // Color avg = sum / (m_output->resolution().x() *
        // m_output->resolution().y()); if(avg.r() == 0.0f)
        //     avg.r() += 0.0001f;
        // if(avg.g() == 0.0f)
        //     avg.g() += 0.0001f;
        // if(avg.b() == 0.0f)
        //     avg.b() += 0.0001f;
        // Color wb = 2.0f * avg;

        Color wb = whitepoint->evaluate(Point2(0, 0), cont);

        for (int y = 0; y < m_output->resolution().y(); y++) {
            for (int x = 0; x < m_output->resolution().x(); x++) {
                Point2i pixel(x, y);
                Color rgb = m_input->get(pixel);

                float maxIntensity = 0.0f;
                bool isEmissive    = false;
                if (rgb.r() > 1.0f || rgb.g() > 1.0f || rgb.b() > 1.0f) {
                    maxIntensity = max(max(rgb.r(), rgb.g()), rgb.b());
                    isEmissive   = true;
                }

                if (isEmissive)
                    rgb /= (maxIntensity + 0.001f);

                // Probably wrong but enables to manipulates channel even if it
                // is 0
                if (rgb.r() == 0.0f)
                    rgb.r() += 0.001f;
                if (rgb.g() == 0.0f)
                    rgb.g() += 0.001f;
                if (rgb.b() == 0.0f)
                    rgb.b() += 0.001f;

                // prevent dividing by 0
                if (wb.r() == 0.0f)
                    wb.r() += 0.001f;
                if (wb.g() == 0.0f)
                    wb.g() += 0.001f;
                if (wb.b() == 0.0f)
                    wb.b() += 0.001f;

                // Color Balancing
                rgb.r() /= wb.r();
                rgb.g() /= wb.g();
                rgb.b() /= wb.b();

                float maxIntensity_WB = 0.0f;
                bool isEmissive_WB    = false;
                if (rgb.r() > 1.0f || rgb.g() > 1.0f || rgb.b() > 1.0f) {
                    maxIntensity_WB = max(max(rgb.r(), rgb.g()), rgb.b());
                    isEmissive_WB   = true;
                }

                if (isEmissive_WB)
                    rgb /= maxIntensity_WB;

                if (isEmissive)
                    rgb *= (maxIntensity + 0.001f);

                m_output->get(pixel) = rgb;
            }
        }

        m_output->save();
    }

    std::string toString() const override {
        return tfm::format(
            "White Balance[\n"
            "  input = %s,\n"
            "  output = %s,\n"
            "  whitepoint = %s, \n"
            "]",
            indent(m_input),
            indent(m_output),
            indent(whitepoint)
        );
    }
};

} // namespace lightwave

REGISTER_POSTPROCESS(Color_Balance, "color_balance");
