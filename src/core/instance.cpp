#include <lightwave/core.hpp>
#include <lightwave/instance.hpp>
#include <lightwave/registry.hpp>
#include <lightwave/sampler.hpp>

namespace lightwave {

void Instance::transformFrame(SurfaceEvent &surf, const Vector &wo,
                              Context &cont, const float time = 0.0f) const {
    SimpleTransform T = m_transform->interpolate(time);
    surf.position     = T.apply(surf.position);

    if (m_normal) {
        const Color normal = 2 * m_normal->evaluate(surf.uv, cont) - Color(1);
        const Vector v     = Vector(normal.r(), normal.g(), normal.b());
        surf.shadingNormal = surf.shadingFrame().toWorld(v);
    }

    surf.shadingNormal = T.applyNormal(surf.shadingNormal).normalized();

    const auto tangent   = surf.tangent.normalized();
    const auto bitangent = surf.geometryNormal.cross(surf.tangent).normalized();
    surf.pdf /= T.apply(tangent).cross(T.apply(bitangent)).length();

    const Vector newGeoNormal = T.applyNormal(surf.geometryNormal);
    surf.tangent              = T.apply(surf.tangent);
    surf.geometryNormal       = newGeoNormal.normalized();
}

inline void validateIntersection(const Intersection &its) {
    // use the following macros to make debugginer easier:
    // * assert_condition(condition, { ... });
    // * assert_normalized(vector, { ... });
    // * assert_ortoghonal(vec1, vec2, { ... });
    // * assert_finite(value or vector or color, { ... });

    // each assert statement takes a block of code to execute when it fails
    // (useful for printing out variables to narrow done what failed)

    assert_finite(its.t, {
        logger(
            EError,
            "  your intersection produced a non-finite intersection distance");
        logger(EError, "  offending shape: %s", its.instance->shape());
    });
    assert_condition(its.t >= Epsilon, {
        logger(EError,
               "  your intersection is susceptible to self-intersections");
        logger(EError, "  offending shape: %s", its.instance->shape());
        logger(EError,
               "  returned t: %.3g (smaller than Epsilon = %.3g)",
               its.t,
               Epsilon);
    });
}

bool Instance::intersect(const Ray &worldRay, Intersection &its, Sampler &rng,
                         Context &cont) const {
    const auto previousIts(its);
    const auto prevAlphaMask = its.alphaMask;
    its.alphaMask            = m_alpha.get();
    if (!m_transform) {
        // fast path, if no transform is needed
        const Ray localRay = worldRay;
        const bool wasIntersected =
            m_shape->intersect(localRay, its, rng, cont);
        if (wasIntersected) {
            its.instance = this;
            validateIntersection(its);
        }
        its.alphaMask = prevAlphaMask;
        return wasIntersected;
    }

    const float previousT = its.t;
    SimpleTransform T     = m_transform->interpolate(worldRay.time);
    Ray localRay          = T.inverse(worldRay);

    const float dLength = localRay.direction.length();
    if (dLength == 0)
        return false;
    localRay.direction /= dLength;

    its.t *= dLength;
    // hints:
    // * transform the ray (do not forget to normalize!)
    // * how does its.t need to change?

    const bool wasIntersected = m_shape->intersect(localRay, its, rng, cont);
    if (!wasIntersected) {
        its.t = previousT;
        its.alphaMask = prevAlphaMask;
        return false;
    }

    if(m_alpha && rng.next() > m_alpha->scalar(its.uv, cont)){
        Ray newR(its.position, localRay.direction);
        float new_t = its.t;

        while(true){
            newR.origin = its.position;
            if(!m_shape->intersect(newR, its, rng, cont) || this != its.instance){
                its = previousIts;
                return  false;
            }
            new_t += its.t;

            if(previousIts.t < new_t){
                its = previousIts;
                return false;
            }

            if(rng.next() <= m_alpha->scalar(its.uv, cont)){
                break;
            }
        }
        its.t = new_t;
    }

    its.instance = this;
    validateIntersection(its);
    assert_finite(its.position, {
        logger(EError,
               "non-finite position, offending shape: %s",
                m_shape->toString());
                logger(EError, "  returned its.t: %g", its.t);
        logger(EError,
                "  for input ray %s in dir %s",
                worldRay.origin,
                worldRay.direction);
    });
    its.t /= dLength;

    transformFrame(its, -localRay.direction, cont, worldRay.time);

    its.alphaMask = prevAlphaMask;
    return true;
}

Bounds Instance::getBoundingBox() const {
    if (!m_transform) {
        // fast path
        return m_shape->getBoundingBox();
    }

    return m_transform->getBoundingBox(m_shape->getBoundingBox());
}

Point Instance::getCentroid() const {
    if (!m_transform) {
        // fast path
        return m_shape->getCentroid();
    }
    SimpleTransform T = m_transform->interpolate(0);
    return T.apply(m_shape->getCentroid());
}

AreaSample Instance::sampleArea(Sampler &rng, Context &cont) const {
    AreaSample sample = m_shape->sampleArea(rng, cont);
    transformFrame(sample, Vector(), cont);
    return sample;
}

} // namespace lightwave

REGISTER_CLASS(Instance, "instance", "default")
