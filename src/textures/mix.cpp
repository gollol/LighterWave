#include <lightwave.hpp>

namespace lightwave {

class MixTexture : public Texture {
    enum class Data_Type { Color, Float, Vector };

    std::vector<std::pair<std::string, Data_Type>> Data_Type_MAPPING = {
        { "RGBA", Data_Type::Color },
        { "FLOAT", Data_Type::Float },
        { "VECTOR", Data_Type::Vector }
    };

    enum class Factor_Vector { Uniform, NonUniform };

    std::vector<std::pair<std::string, Factor_Vector>> Factor_MAPPING = {
        { "UNIFORM", Factor_Vector::Uniform },
        { "NON_UNIFORM", Factor_Vector::NonUniform }
    };

    enum class Blend_Type { Mix, Multiply, Add };

    std::vector<std::pair<std::string, Blend_Type>> Blend_MAPPING = {
        { "MIX", Blend_Type::Mix },
        { "MULTIPLY", Blend_Type::Multiply },
        { "ADD", Blend_Type::Add }
    };

    bool clamp_result_color;
    bool clamp_factor;
    ref<Texture> factor_float;
    ref<Texture> factor_vector;
    ref<Texture> a_color;
    ref<Texture> b_color;
    ref<Texture> a_float;
    ref<Texture> b_float;
    ref<Texture> a_vector;
    ref<Texture> b_vector;

    Data_Type m_type;
    Factor_Vector m_factor_vector;
    Blend_Type m_blend_type;

    Color float_mix(const Point2 &uv, const Context &cont, float fac) const {
        const float a = this->a_float->evaluate(uv, cont).mean();
        const float b = this->b_float->evaluate(uv, cont).mean();

        return Color(b * fac + a * (1.0f - fac));
    }

    // Uniform and Non Uniform is already incorporated into fac
    Color vector_mix(const Point2 &uv, const Context &cont, Color fac) const {
        const Color a = this->a_vector->evaluate(uv, cont);
        const Color b = this->b_vector->evaluate(uv, cont);

        return Color(b * fac + a * (Color(1.0f) - fac));
    }

    Color clamp_color_to_zero_one(Color value) const {
        return Color(clamp(value.r(), 0.0f, 1.0f),
                     clamp(value.g(), 0.0f, 1.0f),
                     clamp(value.b(), 0.0f, 1.0f));
    }

    Color color_mix(const Point2 &uv, const Context &cont, float fac) const {
        const Color a = this->a_color->evaluate(uv, cont);
        const Color b = this->b_color->evaluate(uv, cont);

        Color result;

        switch (m_blend_type) {
        case Blend_Type::Mix:
            result = b * fac + a * (1.0f - fac);
            break;
        case Blend_Type::Multiply:
            result = interpolateLinear(fac, a, a * b);
            break;
        case Blend_Type::Add:
            result = interpolateLinear(fac, a, a + b);
            break;
        default:
            break;
        }

        if (clamp_result_color) {
            result = clamp_color_to_zero_one(result);
        }

        return result;
    }

public:
    MixTexture(const Properties &properties) {
        this->factor_float  = properties.get<Texture>("factor_Float");
        this->factor_vector = properties.get<Texture>("factor_Vector");
        this->a_color       = properties.get<Texture>("a_Color");
        this->b_color       = properties.get<Texture>("b_Color");
        this->a_float       = properties.get<Texture>("a_Float");
        this->b_float       = properties.get<Texture>("b_Float");
        this->a_vector      = properties.get<Texture>("a_Vector");
        this->b_vector      = properties.get<Texture>("b_Vector");

        this->m_type =
            properties.getEnum<Data_Type>("data_type", Data_Type_MAPPING);
        this->m_factor_vector =
            properties.getEnum<Factor_Vector>("factor_mode", Factor_MAPPING);
        this->m_blend_type =
            properties.getEnum<Blend_Type>("blend_type", Blend_MAPPING);

        this->clamp_result_color = properties.get<bool>("clamp_result");
        this->clamp_factor       = properties.get<bool>("clamp_factor");
    }

    /// @note See
    /// https://github.com/iRath96/raymond/blob/fd6341d2c1c63698705343c3aabfffb1c13905a5/raymond/device/nodes/nodes.hpp#L1268
    Color evaluate(const Point2 &uv, const Context &cont) const override {
        switch (m_type) {
        case Data_Type::Float: {
            float fac_float = this->factor_float->evaluate(uv, cont).mean();
            if (clamp_factor) {
                fac_float = clamp(fac_float, 0.0f, 1.0f);
            }
            return float_mix(uv, cont, fac_float);
        }
        case Data_Type::Vector: {
            Color fac_vector;
            if (m_factor_vector == Factor_Vector::Uniform) {
                float mean_vector =
                    this->factor_float->evaluate(uv, cont).mean();
                fac_vector = Color(mean_vector);
            } else {
                fac_vector = this->factor_vector->evaluate(uv, cont);
            }

            if (clamp_factor) {
                fac_vector = clamp_color_to_zero_one(fac_vector);
            }

            return vector_mix(uv, cont, fac_vector);
        }
        case Data_Type::Color: {
            float fac_color = this->factor_float->evaluate(uv, cont).mean();
            if (clamp_factor) {
                fac_color = clamp(fac_color, 0.0f, 1.0f);
            }

            return color_mix(uv, cont, fac_color);
        }
        }

        assert(!"Unsupported Data_Type reached!");
        return Color::black();
    }

    std::string toString() const override {
        return tfm::format(
            "MathTexture[\n"
            "  clamp_factor = %s\n"
            "  clamp_result = %s\n"
            "  data_type = %s\n"
            "  factor_mode = %s\n"
            "  blend_type = %s\n"
            "  factor_float = %s\n"
            "  factor_vector = %s\n"
            "  a_color = %s\n"
            "  b_color = %s\n"
            "  a_float = %s\n"
            "  b_float = %s\n"
            "  a_vector = %s\n"
            "  b_vector = %s\n"
            "]",
            std::to_string(this->clamp_factor),
            std::to_string(this->clamp_result_color),
            enumToString<Data_Type>(this->m_type, Data_Type_MAPPING),
            enumToString<Factor_Vector>(this->m_factor_vector, Factor_MAPPING),
            enumToString<Blend_Type>(this->m_blend_type, Blend_MAPPING),
            indent(this->factor_float),
            indent(this->factor_vector),
            indent(this->a_color),
            indent(this->b_color),
            indent(this->a_float),
            indent(this->b_float),
            indent(this->a_vector),
            indent(this->b_vector));
    }
};

} // namespace lightwave

REGISTER_TEXTURE(MixTexture, "mix")
