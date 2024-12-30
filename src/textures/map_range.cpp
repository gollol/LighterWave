#include <lightwave.hpp>

namespace lightwave {

class MapRangeTexture : public Texture {
    enum class InterpolationMethod {
        Linear,
        SteppedLinear,
        SmoothStep,
        SmootherStep
    };
    std::vector<std::pair<std::string, InterpolationMethod>>
        InterpolationMapping = {
            { "LINEAR", InterpolationMethod::Linear },
            { "STEPPED", InterpolationMethod::SteppedLinear },
            { "SMOOTHSTEP", InterpolationMethod::SmoothStep },
            { "SMOOTHERSTEP", InterpolationMethod::SmootherStep }
        };

    enum class DataType { Float, Vector };
    std::vector<std::pair<std::string, DataType>> TypeMapping = {
        { "FLOAT", DataType::Float }, { "VECTOR", DataType::Vector }
    };

    ref<Texture> m_value;
    ref<Texture> m_fromMin;
    ref<Texture> m_fromMax;
    ref<Texture> m_toMin;
    ref<Texture> m_toMax;

    InterpolationMethod m_interpoationMethod;
    DataType m_type;
    bool do_clamp;

public:
    MapRangeTexture(const Properties &properties) {
        this->m_value   = properties.get<Texture>("value");
        this->m_fromMin = properties.get<Texture>("fromMin");
        this->m_fromMax = properties.get<Texture>("fromMax");
        this->m_toMin   = properties.get<Texture>("toMin");
        this->m_toMax   = properties.get<Texture>("toMax");
        this->do_clamp  = properties.get<bool>("clamp");

        this->m_interpoationMethod = properties.getEnum<InterpolationMethod>(
            "interpolation_type", InterpolationMapping);
        this->m_type = properties.getEnum<DataType>("data_type", TypeMapping);
        logger(EInfo, "map_range initialized");
    }

    /// @note See
    /// https://github.com/iRath96/raymond/blob/fd6341d2c1c63698705343c3aabfffb1c13905a5/raymond/device/nodes/nodes.hpp#L1268
    Color evaluate(const Point2 &uv, const Context &cont) const override {
        const float value   = this->m_value->evaluate(uv, cont).mean();
        const float fromMin = this->m_fromMin->evaluate(uv, cont).mean();
        const float fromMax = this->m_fromMax->evaluate(uv, cont).mean();
        const float toMin   = this->m_toMin->evaluate(uv, cont).mean();
        const float toMax   = this->m_toMax->evaluate(uv, cont).mean();

        // if (rand() % 1000 == 0)
        //     logger(EInfo, "%d", value);

        float out = value;
        switch (m_interpoationMethod) {
        case InterpolationMethod::Linear: {
            out -= fromMin;
            out *= (toMax - toMin) / (fromMax - fromMin);
            out += toMin;
        } break;
        case InterpolationMethod::SteppedLinear:
            break;
        case InterpolationMethod::SmoothStep:
            break;
        case InterpolationMethod::SmootherStep:
            break;
        default:
            break;
        }

        if (do_clamp) {
            out = clamp(out, toMin, toMax);
        }

        return Color(out);
    }

    std::string toString() const override {
        return tfm::format(
            "MathTexture[\n"
            "  interpolation_method = %s\n"
            "  data_type = %s\n"
            "  clamp = %s\n"
            "  value = %s\n"
            "  fromMin = %s\n"
            "  fromMax = %s\n"
            "  toMin = %s\n"
            "  toMax = %s\n"
            "]",
            enumToString<InterpolationMethod>(this->m_interpoationMethod,
                                              InterpolationMapping),
            enumToString<DataType>(this->m_type, TypeMapping),
            std::to_string(this->do_clamp),
            indent(this->m_value),
            indent(this->m_fromMin),
            indent(this->m_fromMax),
            indent(this->m_toMin),
            indent(this->m_toMax));
    }
};

} // namespace lightwave

REGISTER_TEXTURE(MapRangeTexture, "map_range")
