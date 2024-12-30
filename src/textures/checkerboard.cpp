#include <lightwave.hpp>

namespace lightwave {

class Checkerboard : public Texture {
    ref<Texture> m_scale;
    ref<Texture> m_color0;
    ref<Texture> m_color1;

public:
    Checkerboard(const Properties &properties) {
        this->m_scale  = properties.get<Texture>("scale");
        this->m_color0 = properties.get<Texture>("color0");
        this->m_color1 = properties.get<Texture>("color1");
    }

    /// @note See
    /// https://github.com/iRath96/raymond/blob/fd6341d2c1c63698705343c3aabfffb1c13905a5/raymond/device/nodes/nodes.hpp#L247
    Color evaluate(const Point2 &uv, const Context &cont) const override {
        Color scale_color = this->m_scale->evaluate(uv, cont);
        Vector2 scale(scale_color.r(), scale_color.g());
        const auto p =
            Vector2(Point2(0.000001f) + (Vector2(uv) * scale)) * 0.999999f;
        Vector2i idx = Vector2i(floor(p.x()), floor(p.y()));

        const bool which = (idx.x() ^ idx.y()) & 1;
        return which ? this->m_color0->evaluate(uv, cont)
                     : this->m_color1->evaluate(uv, cont);
    }

    std::string toString() const override {
        return tfm::format(
            "Checkerboard[\n"
            "  scale = %s\n"
            "  color0 = %s\n"
            "  color1 = %s\n"
            "]",
            this->m_scale,
            this->m_color0,
            this->m_color1);
    }
};

} // namespace lightwave

REGISTER_TEXTURE(Checkerboard, "checkerboard")
