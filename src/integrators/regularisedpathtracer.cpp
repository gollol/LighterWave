#include <lightwave.hpp>
#include <lightwave/opsr_data.hpp>
#include <lightwave/opsr_path.hpp>

namespace lightwave {

class RegularisedPathTracer : public SamplingIntegrator {
    int m_maxDepth;
    bool m_nee;
    bool m_mis;
    float m_roughnessRes;
    AttenuationFactors m_attFactors;

    float miWeight(float a, float b) const {
        if (std::isinf(a))
            return 1;
        return a / (a + b);
    }

    Color regularisedNextEventEstimation(const Intersection &its,
                                         const float time, Sampler &rng,
                                         float regularisedRoughness,
                                         std::vector<float> &attenuation,
                                         Context &cont) const {
        PROFILE("RNEE")

        auto lightSample = m_scene->sampleLight(rng);
        if (!lightSample)
            return Color(0);

        auto directSample = lightSample.light->sampleDirect(its.position, rng);
        if (m_scene->intersect(Ray(its.position, directSample.wi, time, 0),
                               directSample.distance,
                               rng,
                               cont))
            // occluded
            return Color(0);

        auto B = its.evaluateBsdf(directSample.wi, cont, regularisedRoughness);
        const float misWeight =
            m_mis && lightSample.light->canBeIntersected()
                ? miWeight(lightSample.probability * directSample.pdf, B.pdf)
                : 1;
        return B.value * misWeight * directSample.weight /
               lightSample.probability;
    }

public:
    RegularisedPathTracer(const Properties &properties)
        : SamplingIntegrator(properties) {
        m_maxDepth     = properties.get<int>("depth", 2);
        m_nee          = properties.get<bool>("nee", true);
        m_mis          = properties.get<bool>("mis", false);
        m_roughnessRes = 4;
        auto rougheningMode =
            properties.get<std::string>("rougheningMode", "low");

        if (rougheningMode == "low") {
            m_attFactors = attenuationLow;
        } else if (rougheningMode == "moderate") {
            m_attFactors = attenuationModerate;
        } else if (rougheningMode == "strong") {
            m_attFactors = attenuationStrong;
        } else if (rougheningMode == "aggressive") {
            m_attFactors = attenuationAggressive;
        } else {
            logger(EError, "Invalid roughening mode");
        }
    }

    Color Li(const Ray &cameraRay, Sampler &rng, Context &cont) override {
        PROFILE("Li")

        auto L          = Color(0);
        auto throughput = Color(1);
        auto ray        = cameraRay;

        std::vector<float> pathRoughness;
        std::vector<float> attenuationFactors;
        Ray regularisedRay(cameraRay);

        for (int depth = 1;; depth++) {

            auto its = m_scene->intersect(ray, rng, cont);
            // Direct illumination
            if (depth == 1) {
                if (auto E = its.evaluateEmission(cont)) {
                    L += E.value;
                }
            }

            if (!its || (depth >= m_maxDepth))
                break;

            pathRoughness.push_back(its.vertexRoughness(cont) *
                                    its.vertexRoughness(cont));
            float regularisedRoughness = std::sqrt(
                getRougheningOPSR(m_roughnessRes, pathRoughness, m_attFactors));

            // Regularised NEE
            if (m_nee)
                L += throughput *
                     regularisedNextEventEstimation(its,
                                                    ray.time,
                                                    rng,
                                                    regularisedRoughness,
                                                    attenuationFactors,
                                                    cont);

            auto regularisedBsdfSample =
                its.sampleBsdf(rng, cont, regularisedRoughness);
            if (regularisedBsdfSample) {
                // Shooting a regularised ray to enable MIS
                regularisedRay           = Ray(ray);
                regularisedRay.origin    = its.position;
                regularisedRay.direction = regularisedBsdfSample.wi;

                auto regIts = m_scene->intersect(regularisedRay, rng, cont);
                if (auto E = regIts.evaluateEmission(cont)) {
                    const auto misWeight =
                        !m_nee || (E.pdf == 0) ? 1
                        : m_mis ? miWeight(regularisedBsdfSample.pdf, E.pdf)
                                : 0;
                    L += throughput * misWeight * E.value *
                         regularisedBsdfSample.weight;
                }
            }

            auto B = its.sampleBsdf(rng, cont);
            if (!B)
                break;

            // Use the unbiased path prefix to introduce less further bias
            throughput *= B.weight;
            ray.origin    = its.position;
            ray.direction = B.wi;
        }

        return L;
    }

    std::string toString() const override {
        return tfm::format(
            "RegularisedPathTracer[\n"
            "  sampler = %s,\n"
            "  image = %s,\n"
            "]",
            indent(m_sampler),
            indent(m_image));
    }
};

} // namespace lightwave

REGISTER_INTEGRATOR(RegularisedPathTracer, "regularisedpathtracer")
