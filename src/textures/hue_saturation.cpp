#include <lightwave.hpp>

namespace lightwave {

class HueSaturationTexture : public Texture {
    // Original color value
    ref<Texture> m_original;

    // Modifiers for HSV
    ref<Texture> m_hue;
    ref<Texture> m_saturation;
    ref<Texture> m_value;
    ref<Texture> m_factor;

public:
    HueSaturationTexture(const Properties &properties) {
        this->m_original = properties.getChild<Texture>(true);
        
        this->m_hue = properties.get<Texture>("hue");
        this->m_saturation = properties.get<Texture>("saturation");
        this->m_value = properties.get<Texture>("value");
        this->m_factor = properties.get<Texture>("factor");
    }

    /// @note See https://github.com/iRath96/raymond/blob/fd6341d2c1c63698705343c3aabfffb1c13905a5/raymond/device/nodes/nodes.hpp#L1268 
    Color evaluate(const Point2 &uv, const Context &cont) const override {
        const Color origColor = this->m_original->evaluate(uv, cont);
        Color hsvColor = convertRGBToHSV(origColor);
        
        hsvColor.r() = static_cast<float>(fmod(hsvColor.r() + this->m_hue->scalar(uv, cont) + 0.5f, 1));
        hsvColor.g() = saturate(hsvColor.g() * this->m_saturation->scalar(uv, cont));
        hsvColor.b() *= this->m_value->scalar(uv, cont);

        Color result = max(convertHSVToRGB(hsvColor), Color::black());

        // Interpolate between original color and HSV modified one. "How much the modification should be applied"
        return lerp(origColor, result, this->m_factor->scalar(uv, cont));
    }

    std::string toString() const override {
        return tfm::format("HueSaturationTexture[\n"
                           "  original = %s\n"
                           "  hue = %s\n"
                           "  saturation = %s\n"
                           "  value = %s\n"
                           "  factor = %s\n"
                           "]",
                           indent(this->m_original),
                           indent(this->m_hue),
                           indent(this->m_saturation),
                           indent(this->m_value),
                           indent(this->m_factor));
    }
};

} // namespace lightwave

REGISTER_TEXTURE(HueSaturationTexture, "hue_saturation")
