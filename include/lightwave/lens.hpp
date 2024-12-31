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
    float m_thickness;
    /// @brief Refractive index of the lens.
    float m_eta;
    /// @brief Radius of the lens.
    float m_apertureRadius;
    /// @brief Radius of curvature of the lens.
    float m_curvatureRadius;

public:
    Lens(const Properties &properties) {
        m_thickness       = 0.001f * properties.get<float>("thickness");
        m_eta             = properties.get<float>("eta");
        m_apertureRadius  = 0.001f * properties.get<float>("apertureRadius");
        m_curvatureRadius = 0.001f * properties.get<float>("curvatureRadius");
    }

    /// @brief Calculates the intersection between the given ray and the
    /// spherical lens.
    virtual bool intersect(const Ray &ray, const float zPosition, float &t,
                           Vector &n) const = 0;

    // Getter methods
    float getThickness() const { return m_thickness; }
    float getEta() const { return m_eta; }
    float getApertureRadius() const { return m_apertureRadius; }
    float getCurvatureRadius() const { return m_curvatureRadius; }

    // Setter methods
    void setThickness(const float thickness) { m_thickness = thickness; }
};

} // namespace lightwave
