#include <lightwave.hpp>

namespace lightwave {

class DirectLightIntegrator : public SamplingIntegrator {
public:
    DirectLightIntegrator(const Properties &properties)
        : SamplingIntegrator(properties) {}

    Color Li(const Ray &ray, Sampler &rng, Context &cont) override {
        Intersection its = m_scene->intersect(ray, rng, cont);
        auto L           = Color(0);
        L += its.evaluateEmission(cont).value;
        if (!its)
            return L;

        if (const auto lightSample = m_scene->sampleLight(rng)) {
            if (!lightSample.light->canBeIntersected()) {
                const auto directSample =
                    lightSample.light->sampleDirect(its.position, rng);

                Ray shadowRay       = { ray };
                shadowRay.origin    = its.position;
                shadowRay.direction = directSample.wi;
                if (!m_scene->intersect(shadowRay,
                                        directSample.distance,
                                        rng,
                                        cont)) { // unoccluded
                    auto B = its.evaluateBsdf(shadowRay.direction, cont);
                    L +=
                        B.value * directSample.weight / lightSample.probability;
                }
            }
        }

        if (auto B = its.sampleBsdf(rng, cont)) {
            Ray next       = { ray };
            next.origin    = its.position;
            next.direction = B.wi;
            auto its2      = m_scene->intersect(next, rng, cont);
            L += B.weight * its2.evaluateEmission(cont).value;
        } else {
            // bsdf sampling failed, L is unchanged
        }

        return L;
    }

    std::string toString() const override {
        return tfm::format(
            "DirectLightIntegrator[\n"
            "  sampler = %s,\n"
            "  image = %s,\n"
            "]",
            indent(m_sampler),
            indent(m_image));
    }
};

} // namespace lightwave

REGISTER_INTEGRATOR(DirectLightIntegrator, "direct")
