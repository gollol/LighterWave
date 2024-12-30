#include <lightwave.hpp>

namespace lightwave {

#define SPLIT_COMPONENT_MAPPING {   \
        {"r", SplitComponent::R},   \
        {"g",  SplitComponent::G }, \
        {"b",  SplitComponent::B }, \
    }

#define CONVERT_MODE_MAPPING {          \
        {"NONE", ConvertMode::NONE},    \
        {"HSV",  ConvertMode::HSV },    \
        /* {"HSL",  ConvertMode::B }, */\
    }

class SeparateTexture : public Texture {

    enum class SplitComponent {
        R,
        G,
        B,
    };

    enum class ConvertMode {
        NONE,
        HSV,
        //HSL,  // Currently not supported
    };

    // Original color value
    ref<Texture> m_original;

    // The component to separate
    SplitComponent m_component;

    // What to convert the color into before separating component
    ConvertMode m_convert;
public:
    SeparateTexture(const Properties &properties) {
        this->m_original = properties.getChild<Texture>(true);
        
        this->m_component  = properties.getEnum<SplitComponent>("component", SPLIT_COMPONENT_MAPPING);
        this->m_convert  = properties.getEnum<ConvertMode>("convert", CONVERT_MODE_MAPPING);
    }

    Color evaluate(const Point2 &uv, const Context &cont) const override {
        Color color = this->m_original->evaluate(uv, cont);
        
        // First, check if conversion is needed
        switch (this->m_convert) {
        case ConvertMode::HSV :
            color = convertRGBToHSV(color);
            break;
        default:
            // Do nothing, we're not converting
            break;
        }

        // Now, return separated component as new color
        switch (this->m_component) {
        case SplitComponent::R:
            return Color(color.r());
        case SplitComponent::G:
            return Color(color.g());
        case SplitComponent::B:
            return Color(color.b());
        default:
            // This is here to stop compiler warnings, because the enum should always be set
            return Color(0);
        }
    }

    std::string toString() const override {
        return tfm::format("SeparateTexture[\n"
                           "  original = %s\n"
                           "  component = %s\n"
                           "  convert = %s\n"
                           "]",
                           indent(this->m_original),
                           enumToString<SplitComponent>(m_component, SPLIT_COMPONENT_MAPPING),
                           enumToString<ConvertMode>(m_convert, CONVERT_MODE_MAPPING));
    }
};

} // namespace lightwave

REGISTER_TEXTURE(SeparateTexture, "separate")
