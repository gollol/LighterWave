#include <lightwave.hpp>

namespace lightwave {

class ColorInverting : public Postprocess {

public:
    ColorInverting(const Properties &properties) : Postprocess(properties) 
    {
    }

    /// @brief Colors of input image get inverted. Colors above 1.0f get divided by the highest value r, g, or b - then inverted and again multiplied by the highest value.
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

                //check wether Color is above 1.0f
                float maxIntensity = 0.0f;
                bool isEmissive = false;
                if(rgb.r() > 1.0f || rgb.g() > 1.0f || rgb.b() > 1.0f)
                {
                    maxIntensity = max(max(rgb.r(), rgb.g()), rgb.b());
                    isEmissive = true;
                }

                //inverting
                if(isEmissive)
                {
                    rgb /= maxIntensity;
                    rgb = Color(1) - rgb;
                    rgb *= maxIntensity;
                }
                else
                {
                    rgb = Color(1) - rgb;
                }

                m_output->get(pixel) = rgb;
            }
        }

        m_output->save();        
    }

    std::string toString() const override { return "Color inverting BAMM!"; }
};

} // namespace lightwave

REGISTER_POSTPROCESS(ColorInverting, "color_inverting");
