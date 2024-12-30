#include <lightwave.hpp>

namespace lightwave {

/**
 * @brief A perspective camera with a given field of view angle and transform.
 *
 * In local coordinates (before applying m_transform), the camera looks in
 * positive z direction [0,0,1]. Pixels on the left side of the image ( @code
 * normalized.x < 0 @endcode ) are directed in negative x direction ( @code
 * ray.direction.x < 0 ), and pixels at the bottom of the image ( @code
 * normalized.y < 0 @endcode ) are directed in negative y direction ( @code
 * ray.direction.y < 0 ).
 */
class Perspective : public Camera {
    Vector2 m_span;

public:
    Perspective(const Properties &properties) : Camera(properties) {
        const float fov           = properties.get<float>("fov");
        const std::string fovAxis = properties.get<std::string>("fovAxis");

        const float fovTan = tan(fov * Pi / 360);
        const float aspect = m_resolution.x() / float(m_resolution.y());

        if (fovAxis == "x") {
            m_span = { fovTan, fovTan / aspect };
        } else {
            m_span = { fovTan * aspect, fovTan };
        }
        // hints:
        // * precompute any expensive operations here (most importantly
        // trigonometric functions)
        // * use m_resolution to find the aspect ratio of the image
    }

    CameraSample sample(const Point2 &normalized, Sampler &rng) const override {
        Point origin{ 0, 0, 0 };
        Vector direction{ m_span * Vector2(normalized), 1 };
        const float time  = rng.next();
        SimpleTransform T = m_transform->interpolate(time);

        return {
            .ray    = T.apply(Ray(origin, direction, time, 0)).normalized(),
            .weight = Color::white(),
        };
        // hints:
        // * use m_transform to transform the local camera coordinate system
        // into the world coordinate system
    }

    std::string toString() const override {
        return tfm::format(
            "Perspective[\n"
            "  width = %d,\n"
            "  height = %d,\n"
            "  transform = %s,\n"
            "]",
            m_resolution.x(),
            m_resolution.y(),
            indent(m_transform));
    }
};

} // namespace lightwave

REGISTER_CAMERA(Perspective, "perspective")
