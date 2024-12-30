#include <lightwave.hpp>

namespace lightwave {

#define SAMPLING_MODE_MAPPING                                                  \
    {                                                                          \
        { "BSDF", SamplingMode::BSDF },                                        \
        { "NEE", SamplingMode::NEE },                                          \
        { "MIS", SamplingMode::MIS },                                          \
    }

class PathTracer : public SamplingIntegrator {
    int m_maxDepth;

    enum class SamplingMode {
        BSDF,
        NEE,
        MIS,
    };

    SamplingMode m_samplingMode;
    bool m_fog;
    ref<Texture> m_fogColor;

    float miWeight(float a, float b) const {
        if (std::isinf(a))
            return 1;
        return a / (a + b);
    }

    Color nextEventEstimation(const Intersection &its, const float time,
                              Sampler &rng, Context &cont) const {
        PROFILE("NEE")

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
        auto B = its.evaluateBsdf(directSample.wi, cont);

        const float misWeight =
            (m_samplingMode == SamplingMode::MIS) &&
                    lightSample.light->canBeIntersected()
                ? miWeight(lightSample.probability * directSample.pdf, B.pdf)
                : 1;
        return B.value * misWeight * directSample.weight /
               lightSample.probability;
    }

public:
    PathTracer(const Properties &properties) : SamplingIntegrator(properties) {
        m_maxDepth     = properties.get<int>("depth", 2);
        m_samplingMode = properties.getEnum<SamplingMode>(
            "samplingmode", SamplingMode::NEE, SAMPLING_MODE_MAPPING);
        m_fog      = properties.get<bool>("fog", false);
        m_fogColor = properties.getOptional<Texture>("fogColor");
    }

    Color Li(const Ray &cameraRay, Sampler &rng, Context &cont) override {
        PROFILE("Li")

        auto L           = Color(0);
        auto throughput  = Color(1);
        auto ray         = cameraRay;
        auto prevBsdfPdf = Infinity;

        float distance = 0.0f;

        for (int depth = 1;; depth++) {
            auto its = m_scene->intersect(ray, rng, cont);

            while(its.light() && !its.background)
                its = m_scene->intersect(Ray(its.position, ray.direction), rng, cont);

            if (auto E = its.evaluateEmission(cont)) {
                const auto misWeight =
                    std::isinf(prevBsdfPdf) ||
                            (m_samplingMode == SamplingMode::BSDF) ||
                            (E.pdf == 0)
                        ? 1
                    : (m_samplingMode == SamplingMode::MIS)
                        ? miWeight(prevBsdfPdf, E.pdf)
                        : 0;
                L += throughput * misWeight * E.value;
            }

            if (!its || (depth >= m_maxDepth))
                break;

            distance += its.t;

            if (m_samplingMode != SamplingMode::BSDF)
                L += throughput *
                     nextEventEstimation(its, cameraRay.time, rng, cont);

            auto B = its.sampleBsdf(rng, cont);
            if (!B)
                break;

            throughput *= B.weight;
            prevBsdfPdf   = B.pdf;
            ray           = { ray };
            ray.origin    = its.position;
            ray.direction = B.wi;
        }

        if(!m_fog)
            return L;

        float fog = min(1.f, 3 * sqr(distance / 500));
        // fog = sin(distance / 1000 * Pi2);
        return (1.0f - fog) * L +
            fog * m_fogColor->evaluate(Point2(0), cont);
    }

    std::string toString() const override {
        return tfm::format(
            "PathTracer[\n"
            "  sampler = %s,\n"
            "  image = %s,\n"
            "  fog color = %s,\n"
            "  samplingMode = %s,\n"
            "]",
            indent(m_sampler),
            indent(m_image),
            indent(m_fogColor),
            enumToString<SamplingMode>(m_samplingMode, SAMPLING_MODE_MAPPING));
    }
};

} // namespace lightwave

REGISTER_INTEGRATOR(PathTracer, "pathtracer")
