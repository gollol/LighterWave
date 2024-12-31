#include <lightwave.hpp>

namespace lightwave {

class InvertTexture : public Texture {
    // Original color value
    ref<Texture> m_original;

    // Modifiers for invert
    ref<Texture> m_fac;
public:
    InvertTexture(const Properties &properties) {
        this->m_original = properties.getChild<Texture>(true);
        
        this->m_fac = properties.get<Texture>("fac");
    }

    /// @note See https://github.com/iRath96/raymond/blob/fd6341d2c1c63698705343c3aabfffb1c13905a5/raymond/device/nodes/nodes.hpp#L1268 
    Color evaluate(const Point2 &uv, const Context &cont) const override {
        const Color origColor = this->m_original->evaluate(uv, cont);
        const float fac = this->m_fac->evaluate(uv, cont).mean();
        const Color inverColor = Color{1.0f - origColor.r(), 1.0f - origColor.g(), 1.0f - origColor.b()};
        return inverColor * fac + origColor * (1.f - fac);
    }

    std::string toString() const override {
        return tfm::format("InvertTexture[\n"
                           "  original = %s\n"
                           "  fac = %s\n"
                           "]",
                           indent(this->m_original),
                           indent(this->m_fac));
    }
};

} // namespace lightwave

REGISTER_TEXTURE(InvertTexture, "invert")
