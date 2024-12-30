/**
 * @file warp.hpp
 * @brief Contains functions that map one domain to another.
 */

#pragma once

#include <lightwave/math.hpp>

namespace lightwave {

/**
 * @brief Warps a given point from the unit square ([0,0] to [1,1]) to a unit
 * circle (centered around [0,0] with radius 1), with uniform density given by
 * @code 1 / Pi @endcode .
 * @see Based on
 * http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html
 */
inline Point2 squareToUniformDiskConcentric(const Point2 &sample) {
    float r1 = 2 * sample.x() - 1;
    float r2 = 2 * sample.y() - 1;

    float phi, r;
    if (r1 == 0 && r2 == 0) {
        r   = 0;
        phi = 0;
    } else if (r1 * r1 > r2 * r2) {
        r   = r1;
        phi = Pi4 * (r2 / r1);
    } else {
        r   = r2;
        phi = Pi2 - Pi4 * (r1 / r2);
    }

    float cosPhi = std::cos(phi);
    float sinPhi = std::sin(phi);
    return { r * cosPhi, r * sinPhi };
}

/**
 * @brief Warps a given point from the unit square ([0,0] to [1,1]) to a unit
 * sphere (centered around [0,0,0] with radius 1), with uniform density given by
 * @code 1 / (4 * Pi) @endcode .
 */
inline Vector squareToUniformSphere(const Point2 &sample) {
    float z      = 1 - 2 * sample.y();
    float r      = safe_sqrt(1 - z * z);
    float phi    = 2 * Pi * sample.x();
    float cosPhi = std::cos(phi);
    float sinPhi = std::sin(phi);
    return { r * cosPhi, r * sinPhi, z };
}

/**
 * @brief Warps a given point from the unit square ([0,0] to [1,1]) to a unit
 * hemisphere (centered around [0,0,0] with radius 1, pointing in z direction),
 * with respect to solid angle.
 */
inline Vector squareToUniformHemisphere(const Point2 &sample) {
    Point2 p = squareToUniformDiskConcentric(sample);
    float z  = 1.0f - p.x() * p.x() - p.y() * p.y();
    float s  = sqrt(z + 1.0f);
    return { s * p.x(), s * p.y(), z };
}

/// @brief Returns the density of the @ref squareToUniformHemisphere warping.
inline float uniformHemispherePdf() { return Inv2Pi; }

/**
 * @brief Warps a given point from the unit square ([0,0] to [1,1]) to a unit
 * hemisphere (centered around [0,0,0] with radius 1, pointing in z direction),
 * with density given by @code cos( angle( result, [0,0,1] ) ) @endcode .
 */
inline Vector squareToCosineHemisphere(const Point2 &sample) {
    Point2 p = squareToUniformDiskConcentric(sample);
    float z  = safe_sqrt(1.0f - p.x() * p.x() - p.y() * p.y());
    return { p.x(), p.y(), z };
}

/// @brief Returns the density of the @ref squareToCosineHemisphere warping.
inline float cosineHemispherePdf(const Vector &vector) {
    return InvPi * std::max(vector.z(), float(0));
}

/**
 * @brief Warps a given point from the unit square ([0,0] to [1,1]) to a partial unit
 * hemisphere (centered around [0,0,0] with radius 1, pointing in z direction),
 * with respect to solid angle. The size partial hemisphere is given by the size parameter
 * and is clamped between 0.0-PI/2.0 and determines the max radius of the partial hemisphere.
 */
inline Vector squareToUniformPartialHemisphere(const Point2 &sample, const float max_angle) {
    Point2 p = squareToUniformDiskConcentric(sample);
    float offset = 1.0 - cos(max_angle);
    float offset_sq = sqrt(offset);
    p = Point2(p.x() * offset_sq, p.y() * offset_sq);
    float z  = 1.0f - p.x() * p.x() - p.y() * p.y();
    float s  = sqrt(z + 1.0f);
    return { s * p.x(), s * p.y(), z };
}

/// @brief Returns the density of the @ref squareToUniformPartialHemisphere warping.
inline float uniformPartialHemispherePdf(const Vector &vector, const float max_angle) {
    float cos_angle = vector.dot(Vector(0.0, 0.0, 1.0));
    float angle = acos(cos_angle);
    if (angle > max_angle){
        return 0.0;
    }
    return 1.0 / (-2.0 * Pi * (cos(max_angle) - 1.0)); //Derivation https://www.wolframalpha.com/input?i=integral+from+0+to+2+*+pi+%28integral+from+0+to+a+%28sin%29+dx%29dy
}

/// @brief Returns the density of the @ref squareToUniformPartialHemisphere warping asuming the partial part is hit.
inline float uniformPartialHemispherePdfConst(const float max_angle) {
    return uniformPartialHemispherePdf(Vector(0.0, 0.0, 1.0), max_angle); //Derivation https://www.wolframalpha.com/input?i=integral+from+0+to+2+*+pi+%28integral+from+0+to+a+%28sin%29+dx%29dy
}

inline Vector transformFromZ(const Vector &z_offset, const Vector &normal){
    Vector tangent, bitangent;
    buildOrthonormalBasis(normal, tangent, bitangent);
    return (z_offset.x() * tangent  + z_offset.y() * bitangent + z_offset.z() * normal).normalized();
}

/// @brief given a normalised vector this will return the spherical coordinates
inline Vector2 fromCartesianToSpherical(const Vector &vec){
    float x = vec.x(), y = vec.y(), z = vec.z();
    float inclination = safe_acos(z);
    float azimuth = safe_acos(x / safe_sqrt(x * x + y * y));
    if (y < 0.0){
        azimuth = -azimuth;
    }
    return Vector2(inclination, azimuth);
}


} // namespace lightwave
