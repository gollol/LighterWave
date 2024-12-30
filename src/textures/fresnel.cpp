#include "../bsdfs/fresnel.hpp"
#include <lightwave.hpp>

namespace lightwave {

class FresnelTexture : public Texture {
    // Original color value
    ref<Texture> m_ior;

public:
    FresnelTexture(const Properties &properties) {
        this->m_ior = properties.get<Texture>("ior");
    }

    /// @note See
    /// https://github.com/iRath96/raymond/blob/fd6341d2c1c63698705343c3aabfffb1c13905a5/raymond/device/nodes/nodes.hpp
    Color evaluate(const Point2 &uv, const Context &cont) const override {
        const float ior = this->m_ior->evaluate(uv, cont).mean();
        if (cont.wo.isZero()) {
            return Color(1.0f);
        }
        auto wo         = cont.wo;
        float cosI      = wo.dot(Vector(0.0f, 0.0f, 1.0f));
        bool backfacing = cosI < 0;
        float eta       = max(ior, 1e-5);
        eta             = backfacing ? 1 / eta : eta;

        float fac = fresnelDielectric(cosI, eta);

        return Color(fac);
    }

    std::string toString() const override {
        return tfm::format(
            "FresnelTexture[\n"
            "  ior = %s\n"
            "]",
            indent(this->m_ior));
    }
};

} // namespace lightwave

REGISTER_TEXTURE(FresnelTexture, "fresnel")
