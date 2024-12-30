#include <lightwave.hpp>

namespace lightwave {

class CombineXYZTexture : public Texture {
    // Original color value
    ref<Texture> m_x;
    ref<Texture> m_y;
    ref<Texture> m_z;
public:
    CombineXYZTexture(const Properties &properties) {
        this->m_x = properties.get<Texture>("x");
        this->m_y = properties.get<Texture>("y");
        this->m_z = properties.get<Texture>("z");
    }

    /// @note See https://github.com/iRath96/raymond/blob/fd6341d2c1c63698705343c3aabfffb1c13905a5/raymond/device/nodes/nodes.hpp#L1268 
    Color evaluate(const Point2 &uv, const Context &cont) const override {
        const float x = this->m_x->evaluate(uv, cont).mean();
        const float y = this->m_y->evaluate(uv, cont).mean();
        const float z = this->m_z->evaluate(uv, cont).mean();
        return Color(x, y, z);
    }

    std::string toString() const override {
        return tfm::format("InvertTexture[\n"
                           "  x = %s\n"
                           "  y = %s\n"
                           "  z = %s\n"
                           "]",
                           indent(this->m_x),
                           indent(this->m_y),
                           indent(this->m_z));
    }
};

} // namespace lightwave

REGISTER_TEXTURE(CombineXYZTexture, "combine_xyz")
