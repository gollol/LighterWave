#include <lightwave/bsdf.hpp>
#include <lightwave/color.hpp>
#include <lightwave/emission.hpp>
#include <lightwave/instance.hpp>
#include <lightwave/light.hpp>
#include <lightwave/math.hpp>
#include <lightwave/profiler.hpp>
#include <lightwave/sampler.hpp>

namespace lightwave {

// based on the MESA implementation of the GLU library
std::optional<Matrix4x4> invert(const Matrix4x4 &m) {
    Matrix4x4 inv;
    // clang-format off
    inv(0, 0) = m(1, 1) * m(2, 2) * m(3, 3) - 
               m(1, 1) * m(2, 3) * m(3, 2) - 
               m(2, 1) * m(1, 2) * m(3, 3) + 
               m(2, 1) * m(1, 3) * m(3, 2) +
               m(3, 1) * m(1, 2) * m(2, 3) - 
               m(3, 1) * m(1, 3) * m(2, 2);

    inv(1, 0) = -m(1, 0) * m(2, 2) * m(3, 3) + 
                m(1, 0) * m(2, 3) * m(3, 2) + 
                m(2, 0) * m(1, 2) * m(3, 3) - 
                m(2, 0) * m(1, 3) * m(3, 2) - 
                m(3, 0) * m(1, 2) * m(2, 3) + 
                m(3, 0) * m(1, 3) * m(2, 2);

    inv(2, 0) = m(1, 0) * m(2, 1) * m(3, 3) - 
               m(1, 0) * m(2, 3) * m(3, 1) - 
               m(2, 0) * m(1, 1) * m(3, 3) + 
               m(2, 0) * m(1, 3) * m(3, 1) + 
               m(3, 0) * m(1, 1) * m(2, 3) - 
               m(3, 0) * m(1, 3) * m(2, 1);

    inv(3, 0) = -m(1, 0) * m(2, 1) * m(3, 2) + 
                 m(1, 0) * m(2, 2) * m(3, 1) +
                 m(2, 0) * m(1, 1) * m(3, 2) - 
                 m(2, 0) * m(1, 2) * m(3, 1) - 
                 m(3, 0) * m(1, 1) * m(2, 2) + 
                 m(3, 0) * m(1, 2) * m(2, 1);
    // clang-format on

    const float det = m(0, 0) * inv(0, 0) + m(0, 1) * inv(1, 0) +
                      m(0, 2) * inv(2, 0) + m(0, 3) * inv(3, 0);
    if (det == 0)
        return {};

    // clang-format off
    inv(0, 1) = -m(0, 1) * m(2, 2) * m(3, 3) + 
                m(0, 1) * m(2, 3) * m(3, 2) + 
                m(2, 1) * m(0, 2) * m(3, 3) - 
                m(2, 1) * m(0, 3) * m(3, 2) - 
                m(3, 1) * m(0, 2) * m(2, 3) + 
                m(3, 1) * m(0, 3) * m(2, 2);

    inv(1, 1) = m(0, 0) * m(2, 2) * m(3, 3) - 
               m(0, 0) * m(2, 3) * m(3, 2) - 
               m(2, 0) * m(0, 2) * m(3, 3) + 
               m(2, 0) * m(0, 3) * m(3, 2) + 
               m(3, 0) * m(0, 2) * m(2, 3) - 
               m(3, 0) * m(0, 3) * m(2, 2);

    inv(2, 1) = -m(0, 0) * m(2, 1) * m(3, 3) + 
                m(0, 0) * m(2, 3) * m(3, 1) + 
                m(2, 0) * m(0, 1) * m(3, 3) - 
                m(2, 0) * m(0, 3) * m(3, 1) - 
                m(3, 0) * m(0, 1) * m(2, 3) + 
                m(3, 0) * m(0, 3) * m(2, 1);

    inv(3, 1) = m(0, 0) * m(2, 1) * m(3, 2) - 
                m(0, 0) * m(2, 2) * m(3, 1) - 
                m(2, 0) * m(0, 1) * m(3, 2) + 
                m(2, 0) * m(0, 2) * m(3, 1) + 
                m(3, 0) * m(0, 1) * m(2, 2) - 
                m(3, 0) * m(0, 2) * m(2, 1);

    inv(0, 2) = m(0, 1) * m(1, 2) * m(3, 3) - 
               m(0, 1) * m(1, 3) * m(3, 2) - 
               m(1, 1) * m(0, 2) * m(3, 3) + 
               m(1, 1) * m(0, 3) * m(3, 2) + 
               m(3, 1) * m(0, 2) * m(1, 3) - 
               m(3, 1) * m(0, 3) * m(1, 2);

    inv(1, 2) = -m(0, 0) * m(1, 2) * m(3, 3) + 
                m(0, 0) * m(1, 3) * m(3, 2) + 
                m(1, 0) * m(0, 2) * m(3, 3) - 
                m(1, 0) * m(0, 3) * m(3, 2) - 
                m(3, 0) * m(0, 2) * m(1, 3) + 
                m(3, 0) * m(0, 3) * m(1, 2);

    inv(2, 2) = m(0, 0) * m(1, 1) * m(3, 3) - 
                m(0, 0) * m(1, 3) * m(3, 1) - 
                m(1, 0) * m(0, 1) * m(3, 3) + 
                m(1, 0) * m(0, 3) * m(3, 1) + 
                m(3, 0) * m(0, 1) * m(1, 3) - 
                m(3, 0) * m(0, 3) * m(1, 1);

    inv(3, 2) = -m(0, 0) * m(1, 1) * m(3, 2) + 
                 m(0, 0) * m(1, 2) * m(3, 1) + 
                 m(1, 0) * m(0, 1) * m(3, 2) - 
                 m(1, 0) * m(0, 2) * m(3, 1) - 
                 m(3, 0) * m(0, 1) * m(1, 2) + 
                 m(3, 0) * m(0, 2) * m(1, 1);

    inv(0, 3) = -m(0, 1) * m(1, 2) * m(2, 3) + 
                m(0, 1) * m(1, 3) * m(2, 2) + 
                m(1, 1) * m(0, 2) * m(2, 3) - 
                m(1, 1) * m(0, 3) * m(2, 2) - 
                m(2, 1) * m(0, 2) * m(1, 3) + 
                m(2, 1) * m(0, 3) * m(1, 2);

    inv(1, 3) = m(0, 0) * m(1, 2) * m(2, 3) - 
               m(0, 0) * m(1, 3) * m(2, 2) - 
               m(1, 0) * m(0, 2) * m(2, 3) + 
               m(1, 0) * m(0, 3) * m(2, 2) + 
               m(2, 0) * m(0, 2) * m(1, 3) - 
               m(2, 0) * m(0, 3) * m(1, 2);

    inv(2, 3) = -m(0, 0) * m(1, 1) * m(2, 3) + 
                 m(0, 0) * m(1, 3) * m(2, 1) + 
                 m(1, 0) * m(0, 1) * m(2, 3) - 
                 m(1, 0) * m(0, 3) * m(2, 1) - 
                 m(2, 0) * m(0, 1) * m(1, 3) + 
                 m(2, 0) * m(0, 3) * m(1, 1);

    inv(3, 3) = m(0, 0) * m(1, 1) * m(2, 2) - 
                m(0, 0) * m(1, 2) * m(2, 1) - 
                m(1, 0) * m(0, 1) * m(2, 2) + 
                m(1, 0) * m(0, 2) * m(2, 1) + 
                m(2, 0) * m(0, 1) * m(1, 2) - 
                m(2, 0) * m(0, 2) * m(1, 1);
    // clang-format on

    return inv * (1.f / det);
}

// https://mo.mathematik.uni-stuttgart.de/inhalt/beispiel/beispiel1113/
std::optional<Matrix3x3> invert(const Matrix3x3 &m) {
    // Based on Sarrus
    const float det =
        m(0, 0) * m(1, 1) * m(2, 2) + m(0, 1) * m(1, 2) * m(2, 0) +
        m(0, 2) * m(1, 0) * m(2, 1) - m(0, 1) * m(1, 0) * m(2, 2) -
        m(0, 2) * m(1, 1) * m(2, 0) - m(0, 0) * m(1, 2) * m(2, 1);
    if (det == 0)
        return {};

    // Based on Cramer's rule
    Matrix3x3 inv;
    inv(0, 0) = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1);
    inv(0, 1) = -m(0, 1) * m(2, 2) + m(0, 2) * m(2, 1);
    inv(0, 2) = m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1);

    inv(1, 0) = -m(1, 0) * m(2, 2) + m(1, 2) * m(2, 0);
    inv(1, 1) = m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0);
    inv(1, 2) = -m(0, 0) * m(1, 2) + m(0, 2) * m(1, 0);

    inv(2, 0) = m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0);
    inv(2, 1) = -m(0, 0) * m(2, 1) + m(0, 1) * m(2, 0);
    inv(2, 2) = m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);

    return 1.f / det * inv;
}

void buildOrthonormalBasis(const Vector &a, Vector &b, Vector &c) {
    if (abs(a.x()) > abs(a.y())) {
        auto invLen = 1 / sqrt(a.x() * a.x() + a.z() * a.z());
        c           = { a.z() * invLen, 0, -a.x() * invLen };
    } else {
        auto invLen = 1 / sqrt(a.y() * a.y() + a.z() * a.z());
        c           = { 0, a.z() * invLen, -a.y() * invLen };
    }
    b = c.cross(a);
}

EmissionEval Intersection::evaluateEmission(Context &cont) const {
    if (!instance) {
        if (background) {
            // nothing was hit, but a background light is available
            auto E = background->evaluate(-wo);
            E.pdf *= lightProbability;
            E.value /= backgroundProbability;
            return E;
        }

        // nothing was hit, no background light available
        return EmissionEval::invalid();
    }

    if (!instance->emission()) {
        // something was hit, but has no emission
        return EmissionEval::invalid();
    }

    // something was hit and has an emission
    auto E =
        instance->emission()->evaluate(uv, shadingFrame().toLocal(wo), cont);
    return {
        .value = E.value,
        .pdf   = lightProbability * directPdf(),
    };
}

float Intersection::directPdf() const {
    return pdf * sqr(t) / Frame::absCosTheta(shadingFrame().toLocal(wo));
}

BsdfSample Intersection::sampleBsdf(Sampler &rng, Context &cont,
                                    float roughness) const {
    PROFILE("Sample Bsdf")

    if (!instance || !instance->bsdf())
        return BsdfSample::invalid();
    assert_normalized(wo, {});
    auto bsdfSample = instance->bsdf()->sample(
        uv, shadingFrame().toLocal(wo), rng, cont, roughness);
    if (bsdfSample.isInvalid())
        return bsdfSample;
    assert_normalized(bsdfSample.wi, {
        logger(EError, "offending BSDF: %s", instance->bsdf()->toString());
        logger(EError, "  input was: %s with length %f", wo, wo.length());
    });
    bsdfSample.wi = shadingFrame().toWorld(bsdfSample.wi);
    assert_normalized(bsdfSample.wi, {});
    return bsdfSample;
}

BsdfEval Intersection::evaluateBsdf(const Vector &wi, Context &cont,
                                    float roughness) const {
    PROFILE("Evaluate Bsdf")

    if (!instance || !instance->bsdf())
        return BsdfEval::invalid();
    auto local_wo = shadingFrame().toLocal(wo);
    cont.wo       = local_wo;
    return instance->bsdf()->evaluate(uv,
                                      shadingFrame().toLocal(wo),
                                      shadingFrame().toLocal(wi),
                                      cont,
                                      roughness);
}

float Intersection::vertexRoughness(const Context &cont) const {
    PROFILE("Getting Vertex Roughness")

    if (!instance || !instance->bsdf())
        return 0;
    return instance->bsdf()->getRoughness(uv, cont);
}

bool Intersection::testAlpha(const Point2 &testUV, Sampler &rng,
                             Context &cont) const {
    if (!alphaMask)
        return true;
    return alphaMask->scalar(testUV, cont) > rng.next();
}

Light *Intersection::light() const {
    if (!instance)
        return background;
    return instance->light();
}

void Frame::validate() const {
    assert_normalized(normal, {});
    assert_normalized(tangent, {});
    assert_normalized(bitangent, {});
    assert_orthogonal(normal, tangent, {});
    assert_orthogonal(normal, bitangent, {});
    assert_orthogonal(tangent, bitangent, {});
}

} // namespace lightwave
