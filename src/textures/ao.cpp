#include <lightwave.hpp>

namespace lightwave {

class AoTexture : public Texture {
    int m_samples;
    float m_distance;
    Color m_color;

public:
    AoTexture(const Properties &properties) {
        m_samples  = properties.get<int>("samples", 16);
        m_distance = properties.get<float>("distance", 1.0f);
        m_color    = properties.get<Color>("color", Color::white());
    }

    Color evaluate(const Point2 &uv, const Context &ctx) const override {
        int unoccluded = 0;

        Context copy = Context(ctx);

        for (int sample = 0; sample < m_samples; sample++) {
            Vector localDirection =
                squareToUniformHemisphere(ctx.rng->next2D());
            Vector worldDirection =
                ctx.its->shadingFrame().toWorld(localDirection);

            // TODO: add the correct time parameter below. This version should
            // work, as long as the object the texture is applied on is not
            // scaled in a non-uniform fashion over time.
            Ray occlusionRay = Ray(ctx.its->position, worldDirection, 0, 0);
            if (!ctx.scene->intersect(
                    occlusionRay, m_distance, *ctx.rng, copy)) {
                unoccluded++;
            }
        }

        float ao = (float) unoccluded / (float) m_samples;

        return m_color * ao;
    }

    virtual float scalar(const Point2 &uv, const Context &ctx) const override {
        int unoccluded = 0;

        Context copy = Context(ctx);

        for (int sample = 0; sample < m_samples; sample++) {
            Vector localDirection =
                squareToUniformHemisphere(ctx.rng->next2D());
            Vector worldDirection =
                ctx.its->shadingFrame().toWorld(localDirection);

            // TODO: add the correct time parameter below. This version should
            // work, as long as the object the texture is applied on is not
            // scaled in a non-uniform fashion over time.
            Ray occlusionRay = Ray(ctx.its->position, worldDirection, 0, 0);
            if (!ctx.scene->intersect(
                    occlusionRay, m_distance, *ctx.rng, copy)) {
                unoccluded++;
            }
        }

        float ao = (float) unoccluded / (float) m_samples;

        return ao;
    }

    std::string toString() const override {
        return tfm::format(
            "AoTexture[\n"
            "]");
    }
};

} // namespace lightwave

REGISTER_TEXTURE(AoTexture, "ao")
