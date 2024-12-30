#include <lightwave.hpp>

namespace lightwave {

class MathTexture : public Texture {
    enum class Operation { Add, Sub, Mul, MulAdd, Min, Max, Power, Greater, Less};

    std::vector<std::pair<std::string, Operation>> Operations_MAPPING = {
        { "ADD", Operation::Add },      { "SUBTRACT", Operation::Sub },
        { "MULTIPLY", Operation::Mul }, { "MULTIPLY_ADD", Operation::MulAdd },
        { "MINIMUM", Operation::Min },  { "MAXIMUM", Operation::Max },
        { "POWER", Operation::Power }, { "GREATER_THAN", Operation::Greater},
        { "LESS_THAN", Operation::Less}
    };

    ref<Texture> m_value;
    ref<Texture> m_value_001;
    ref<Texture> m_value_002;

    bool do_clamp;
    Operation m_operation;

public:
    MathTexture(const Properties &properties) {
        this->m_value     = properties.get<Texture>("value");
        this->m_value_001 = properties.get<Texture>("value_001");
        this->m_value_002 = properties.get<Texture>("value_002");
        this->do_clamp    = properties.get<bool>("clamp");

        this->m_operation =
            properties.getEnum<Operation>("operation", Operations_MAPPING);
    }

    /// @note See
    /// https://github.com/iRath96/raymond/blob/fd6341d2c1c63698705343c3aabfffb1c13905a5/raymond/device/nodes/nodes.hpp#L1268
    Color evaluate(const Point2 &uv, const Context &cont) const override {
        const float v     = this->m_value->evaluate(uv, cont).mean();
        const float v_001 = this->m_value_001->evaluate(uv, cont).mean();
        const float v_002 = this->m_value_002->evaluate(uv, cont).mean();

        float out = 0.0;
        switch (m_operation) {
        case Operation::Add:
            out = v + v_001;
            break;
        case Operation::Sub:
            out = v - v_001;
            break;
        case Operation::Mul:
            out = v * v_001;
            break;
        case Operation::MulAdd:
            out = v * v_001 + v_002;
            break;
        case Operation::Min:
            out = std::min(v, v_001);
            break;
        case Operation::Max:
            out = std::max(v, v_001);
            break;
        case Operation::Power:
            out = safe_powf(v, v_001);
            break;
        case Operation::Less:
            out = v < v_001 ? 1.0f : 0.0f;
            break;
        case Operation::Greater:
            out = v > v_001 ? 1.0f : 0.0f;
            break;
        default:
            break;
        }

        if (do_clamp) {
            out = clamp(out, 0.0, 1.0);
        }

        return Color(out);
    }

    std::string toString() const override {
        return tfm::format(
            "MathTexture[\n"
            "  operation = %s\n"
            "  clamp = %s\n"
            "  value = %s\n"
            "  value_001 = %s\n"
            "  value_002 = %s\n"
            "]",
            enumToString<Operation>(this->m_operation, Operations_MAPPING),
            std::to_string(this->do_clamp),
            indent(this->m_value),
            indent(this->m_value_001),
            indent(this->m_value_002));
    }
};

} // namespace lightwave

REGISTER_TEXTURE(MathTexture, "math")
