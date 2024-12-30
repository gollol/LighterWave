#include <lightwave.hpp>

namespace lightwave {

class BrightnessContrast : public Postprocess {
    float contrast = 1.0f;
    float brightness = 0.0f;

public:
    BrightnessContrast(const Properties &properties) : Postprocess(properties) 
    {
        contrast = properties.get<float>("contrast");
        brightness = properties.get<float>("brightness");
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

                //check wether Color is above 1.0f
                float maxIntensity = 0.0f;
                bool isEmissive = false;
                if(rgb.r() > 1.0f || rgb.g() > 1.0f || rgb.b() > 1.0f)
                {
                    maxIntensity = max(max(rgb.r(), rgb.g()), rgb.b());
                    isEmissive = true;
                }

                //Blender does not clamp
                rgb += Color(brightness);
                rgb.r() = clamp(rgb.r(), 0.0f, Infinity);
                rgb.g() = clamp(rgb.g(), 0.0f, Infinity);
                rgb.b() = clamp(rgb.b(), 0.0f, Infinity);

                //affine color(grey) transform
                rgb -= Color(0.5f);
                if(contrast == 1.0f)
                {
                    if(rgb.r() <= 0.0f)
                        rgb.r() = -0.5f;
                    if(rgb.r() > 0.0f)
                        rgb.r() = 0.5f;

                    if(rgb.g() <= 0.0f)
                        rgb.g() = -0.5f;
                    if(rgb.g() > 0.0f)
                        rgb.g ()= 0.5f;

                    if(rgb.b() <= 0.0f)
                        rgb.b() = -0.5f;
                    if(rgb.b() > 0.0f)
                        rgb.b() = 0.5f;
                }
                else
                {
                    rgb *= (contrast + 1) / (1 - contrast);

                    if(rgb.r() <= -0.5f)
                        rgb.r() = -0.5f;
                    if(rgb.r() > 0.5f)
                        rgb.r() = 0.5f;

                    if(rgb.g() <= -0.5f)
                        rgb.g() = -0.5f;
                    if(rgb.g() > 0.5f)
                        rgb.g ()= 0.5f;

                    if(rgb.b() <= -0.5f)
                        rgb.b() = -0.5f;
                    if(rgb.b() > 0.5f)
                        rgb.b() = 0.5f;
                }
                
                rgb += Color(0.5f);

                m_output->get(pixel) = rgb;
            }
        }

        m_output->save();        
    }

    std::string toString() const override {
        return tfm::format(
            "Brightness / Contrast[\n"
            "  input = %s,\n"
            "  output = %s,\n"
            "  brightness = %s, \n"
            "  contrast = %s, \n"
            "]",
            indent(m_input),
            indent(m_output),
            indent(brightness),
            indent(contrast)
        );
    }
};

} // namespace lightwave

REGISTER_POSTPROCESS(BrightnessContrast, "brightness_contrast");
