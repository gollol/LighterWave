#include <lightwave.hpp>

namespace lightwave {

class HueSaturationValue : public Postprocess {
    float hue = 0.0f;
    float saturation = 1.0f;
    float value = 1.0f;

public:
    HueSaturationValue(const Properties &properties) : Postprocess(properties) 
    {
        hue = properties.get<float>("hue");
        saturation = properties.get<float>("saturation");
        value = properties.get<float>("value");
    }

    /// @brief Given a Color rgb the pointers h,s,v contain the converted color.
    void rgb2hsv(Color rgb, float* h, float* s, float* v)
    {
        float maxIntensity_ = max(max(rgb.r(), rgb.g()), rgb.b());
        float minIntensity_ = min(min(rgb.r(), rgb.g()), rgb.b());

        float delta = maxIntensity_ - minIntensity_;

        unsigned int maxRGB = 0;
        if(maxIntensity_ == rgb.r())
            maxRGB = 1;
        else if(maxIntensity_ == rgb.g())
            maxRGB = 2;
        else if(maxIntensity_ == rgb.b())
            maxRGB = 3;
        
        switch (maxRGB)
        {
        case 1:
            *h = 60 * fmod(((rgb.g() - rgb.b()) / delta), 6.0f);
            break;
        case 2:
            *h = 60 * ((rgb.b() - rgb.r()) / delta + 2.0f);
            break;
        case 3:
            *h = 60 * ((rgb.r() - rgb.g()) / delta + 4.0f);
            break;
        default:
            break;
        }

        if(*h < 0)
            *h = *h + 360.0f;

        if(delta == 0)
            *h = 0;

        *s = maxIntensity_ == 0.0f ? 0.0f : delta / maxIntensity_;
        *v = maxIntensity_;
    }

    /// @brief Given h,s,v the color pointer rgb contains the converted color.
    void hsv2rgb(float h, float s, float v, Color* rgb)
    {
        float c = v * s;
        float x = c * (1.0f - fabs(fmod(h / 60.0f, 2.0f) - 1.0f));
        float m = v - c;

        Color rgb_ = Color(-1);
        if(h >= 0.0f && h < 60.0f)
            rgb_ = Color(c, x, 0);
        else if(h >= 60.0f && h < 120.0f)
            rgb_ = Color(x, c, 0);
        else if(h >= 120.0f && h < 180.0f)
            rgb_ = Color(0, c, x);
        else if(h >= 180.0f && h < 240.0f)
            rgb_ = Color(0, x, c);
        else if(h >= 240.0f && h < 300.0f)
            rgb_ = Color(x, 0, c);
        else if(h >= 300.0f && h < 360.0f)
            rgb_ = Color(c, 0, x);

        *rgb = Color(rgb_.r() + m, rgb_.g() + m, rgb_.b() + m);
    }

    /// @brief Change the Hue, Saturation and Value of the colors in the image.
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

                if(isEmissive)
                    rgb /= maxIntensity;

                float h; 
                float s; 
                float v;
                rgb2hsv(rgb, &h, &s, &v);

                //modify hsv--------------------
                h = h + (hue - 0.5f) * 360.0f;
                if(h < 0)
                    h = h + 360.0f;
                if(h >= 360.0f)
                    h = h - 360.0f;
                
                s = s * saturation;

                v = v * value;
                //modify hsv--------------------

                hsv2rgb(h, s, v, &rgb);       

                if(isEmissive)
                    rgb *= maxIntensity;

                m_output->get(pixel) = rgb;
            }
        }

        m_output->save();        
    }

    std::string toString() const override {
        return tfm::format(
            "HSV[\n"
            "  input = %s,\n"
            "  output = %s,\n"
            "  hue = %s, \n"
            "  saturation = %s, \n"
            "  value = %s, \n"
            "]",
            indent(m_input),
            indent(m_output),
            indent(hue),
            indent(saturation),
            indent(value)
        );
    }
};

} // namespace lightwave

REGISTER_POSTPROCESS(HueSaturationValue, "hue_saturation_value");
