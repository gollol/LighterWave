#include <lightwave.hpp>

namespace lightwave {
/**
 * @brief Implements a thin lens camera model as described the PBRT book:
 * https://pbr-book.org/3ed-2018/Camera_Models/Projective_Camera_Models#TheThinLensModelandDepthofField
 */
class ThinLens : public Camera {

protected:
    /// @brief Factor by which to scale the x coordiante of the image.
    float scaleX;
    /// @brief Factor by which to scale the y coordinate of the image.
    float scaleY;
    /// @brief The axis of the field of view.
    std::string fovAxis;
    /// @brief The aspect ratio of the image.
    float aspect;
    /// @brief The radius of the lens.
    float radius;
    /// @brief The distance from the lens to the plane of focus.
    float focalDistance;

public:
    ThinLens(const Properties &properties) : Camera(properties) {
        const float fov   = properties.get<float>("fov");
        const float s_fov = std::tan((fov * (Pi / 180.f)) * 0.5f);
        aspect = static_cast<float>(m_resolution.x()) / m_resolution.y();

        fovAxis       = properties.get<std::string>("fovAxis");
        radius        = properties.get<float>("radius");
        focalDistance = properties.get<float>("focalDistance");

        if (fovAxis == "x") {
            scaleX = s_fov;
            scaleY = s_fov * 1.f / aspect;
        } else {
            scaleX = s_fov * aspect;
            scaleY = s_fov;
        }
    }

    CameraSample sample(const Point2 &normalized, Sampler &rng) const override {
        // Sample a point on the lens.
        Point2 pLens = squareToUniformDiskConcentric(rng.next2D());
        pLens.x()    = pLens.x() * radius;
        pLens.y()    = pLens.y() * radius;
        // Get the pixel on the image plane.
        float x = normalized.x() * scaleX;
        float y = normalized.y() * scaleY;

        const float time = rng.next();
        
        Ray cameraRay   = Ray(Vector(0, 0, 0), Vector(x, y, 1), time, 0);
        Point rayOrigin = Point(pLens.x(), pLens.y(), 0);
        Point pFocus    = cameraRay(focalDistance);
        cameraRay = Ray(rayOrigin, pFocus - rayOrigin, time, 0);

        SimpleTransform T = m_transform->interpolate(time);
        cameraRay = T.apply(cameraRay);
        
        return CameraSample{ .ray    = cameraRay.normalized(),
                             .weight = Color(1.0f) };
    }

    std::string toString() const override {
        return tfm::format(
            "ThinLens[\n"
            "  width = %d,\n"
            "  height = %d,\n"
            "  fovAxis = %s,\n"
            "  aspect = %f,\n"
            "  radius = %f,\n"
            "  f = %f,\n"
            "  transform = %s,\n"
            "]",
            m_resolution.x(),
            m_resolution.y(),
            fovAxis,
            aspect,
            radius,
            focalDistance,
            indent(m_transform));
    }
};
} // namespace lightwave

REGISTER_CAMERA(ThinLens, "thinlens")
