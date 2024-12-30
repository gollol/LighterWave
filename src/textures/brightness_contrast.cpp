#include <lightwave.hpp>

namespace lightwave {

class BrightnessContrastTexture : public Texture {
    // Original color value
    ref<Texture> m_original;

    // Modifiers for invert
    ref<Texture> m_bright;
    ref<Texture> m_contrast;

public:
    BrightnessContrastTexture(const Properties &properties) {
        this->m_original = properties.getChild<Texture>(true);

        this->m_bright   = properties.get<Texture>("bright");
        this->m_contrast = properties.get<Texture>("contrast");
    }

    /// @note See
    /// https://github.com/iRath96/raymond/blob/fd6341d2c1c63698705343c3aabfffb1c13905a5/raymond/device/nodes/nodes.hpp#L1268
    Color evaluate(const Point2 &uv, const Context &cont) const override {
        const Color origColor = this->m_original->evaluate(uv, cont);
        const float bright    = this->m_bright->evaluate(uv, cont).mean();
        const float contrast  = this->m_contrast->evaluate(uv, cont).mean();

        const float a = 1.0f + contrast;
        const float b = bright - contrast * 0.5f;

        Color out = Color(max(a * origColor.r() + b, 0.0f),
                          max(a * origColor.g() + b, 0.0f),
                          max(a * origColor.b() + b, 0.0f));

        return out;
    }

    std::string toString() const override {
        return tfm::format(
            "BrightnessContrastTexture[\n"
            "  original = %s\n"
            "  bright = %s\n"
            "  contrast = %s\n"
            "]",
            indent(this->m_original),
            indent(this->m_bright),
            indent(this->m_contrast));
    }
};

} // namespace lightwave

REGISTER_TEXTURE(BrightnessContrastTexture, "brightnesscontrast")
