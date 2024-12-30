#include <lightwave.hpp>

namespace lightwave {

class Exposure : public Postprocess {
    float exposure = 0.0f;

public:
    Exposure(const Properties &properties) : Postprocess(properties) 
    {
        exposure = properties.get<float>("exposure");
    }

    /// @brief Image gets more/less contrast/brightness
    void execute() override 
    {
        //Output initialized to resolution of input
        m_output->initialize(m_input->resolution());

        //Loop over each Pixel
        for(int y = 0; y < m_output->resolution().y(); y++)
        {
            for(int x = 0; x < m_output->resolution().x(); x++)
            {
                Point2i pixel(x,y);
                Color rgb = m_input->get(pixel);

                //exposure change
                rgb *= powf(2, exposure);

                m_output->get(pixel) = rgb;
            }
        }

        m_output->save();        
    }

    std::string toString() const override {
        return tfm::format(
            "Exposure[\n"
            "  input = %s,\n"
            "  output = %s,\n"
            "  exposure = %s, \n"
            "]",
            indent(m_input),
            indent(m_output),
            indent(exposure)
        );
    }
};

} // namespace lightwave

REGISTER_POSTPROCESS(Exposure, "exposure");
