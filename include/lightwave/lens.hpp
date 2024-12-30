/**
 * @file lens.hpp
 * @brief Contains the Lens interface and related structures.
 */

#pragma once

#include <lightwave/core.hpp>
#include <lightwave/properties.hpp>

namespace lightwave {

class Lens : public Object {

protected:
    /// @brief Thickness of the lens.
    float thickness;
    /// @brief Refractive index of the lens.
    float eta;
    /// @brief Radius of the lens.
    float apertureRadius;
    /// @brief Radius of curvature of the lens.
    float curvatureRadius;

public:
    Lens(const Properties &properties) {
        thickness       = 0.001 * properties.get<float>("thickness");
        eta             = properties.get<float>("eta");
        apertureRadius  = 0.001 * properties.get<float>("apertureRadius");
        curvatureRadius = 0.001 * properties.get<float>("curvatureRadius");
    }

    /// @brief Calculates the intersection between the given ray and the
    /// spherical lens.
    virtual bool intersect(const Ray &ray, const float zPosition, float &t,
                           Vector &n) const = 0;

    // Getter methods
    float getThickness() const { return thickness; }
    float getEta() const { return eta; }
    float getApertureRadius() const { return apertureRadius; }
    float getCurvatureRadius() const { return curvatureRadius; }

    // Setter methods
    void setThickness(float thickness) { this->thickness = thickness; }
};

} // namespace lightwave
