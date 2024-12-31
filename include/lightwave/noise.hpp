/**
 * @file noise.hpp
 * @brief Contains functions related to generating noise.
 * @note Adapted from
 * https://projects.blender.org/blender/blender/src/commit/4651f8a08faa578ff1a488e68fe8ee659cd31424/intern/cycles/kernel/svm/noise.h
 */

#pragma once

#include <lightwave/math.hpp>

namespace lightwave {

/**
 * Perlin Noise
 */

/// @brief Fade function / ease curve for transitioning between gradients
/// @param t
/// @return Interpolated / eased value
inline float fade(float t) {
    return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
}

/// @brief Negates the value, if the condition is true and leaves it alone if
/// false
inline float negate_if(float val, int condition) {
    return condition ? -val : val;
}

inline float grad1(int hash, float x) {
    int h   = hash & 15;
    float g = 1.f + static_cast<float>(h & 7);
    return negate_if(g, h & 8) * x;
}

inline float perlin_1d(float x) {
    int X;
    float fx = floorfrac(x, X);
    float u  = fade(fx);

    using namespace lightwave::hash;

    return mix(grad1(hash_uint(X), fx), grad1(hash_uint(X + 1), fx * 1.0f), u);
}

inline float grad2(int hash, float x, float y) {
    int h   = hash & 7;
    float u = h < 4 ? x : y;
    float v = 2.0f * (h < 4 ? y : x);
    return negate_if(u, h & 1) + negate_if(v, h & 2);
}

/**
 * Bilinear Interpolation:
 *
 * v2          v3
 *  @ + + + + @       y
 *  +         +       ^
 *  +         +       |
 *  +         +       |
 *  @ + + + + @       @------> x
 * v0          v1
 *
 */
inline float bi_mix(float v0, float v1, float v2, float v3, float x, float y) {
    float x1 = 1.0f - x;
    return (1.0f - y) * (v0 * x1 + v1 * x) + y * (v2 * x1 + v3 * x);
}

inline float perlin_2d(float x, float y) {
    int X;
    int Y;

    float fx = floorfrac(x, X);
    float fy = floorfrac(y, Y);

    float u = fade(fx);
    float v = fade(fy);

    using namespace lightwave::hash;

    return bi_mix(grad2(hash_uint2(X, Y), fx, fy),
                  grad2(hash_uint2(X + 1, Y), fx - 1.0f, fy),
                  grad2(hash_uint2(X, Y + 1), fx, fy - 1.0f),
                  grad2(hash_uint2(X + 1, Y + 1), fx - 1.0f, fy - 1.0f),
                  u,
                  v);
}

inline float grad3(int hash, float x, float y, float z) {
    int h    = hash & 15;
    float u  = h < 8 ? x : y;
    float vt = ((h == 12) || (h == 14)) ? x : z;
    float v  = h < 4 ? y : vt;
    return negate_if(u, h & 1) + negate_if(v, h & 2);
}

/**
 * Trilinear Interpolation:
 *
 *   v6               v7
 *     @ + + + + + + @
 *     +\            +\
 *     + \           + \
 *     +  \          +  \
 *     +   \ v4      +   \ v5
 *     +    @ + + + +++ + @          z
 *     +    +        +    +      y   ^
 *  v2 @ + +++ + + + @ v3 +       \  |
 *      \   +         \   +        \ |
 *       \  +          \  +         \|
 *        \ +           \ +          +---------> x
 *         \+            \+
 *          @ + + + + + + @
 *        v0               v1
 */
inline float tri_mix(float v0, float v1, float v2, float v3, float v4, float v5,
                     float v6, float v7, float x, float y, float z) {
    float x1 = 1.0f - x;
    float y1 = 1.0f - y;
    float z1 = 1.0f - z;
    return z1 * (y1 * (v0 * x1 + v1 * x) + y * (v2 * x1 + v3 * x)) +
           z * (y1 * (v4 * x1 + v5 * x) + y * (v6 * x1 + v7 * x));
}

inline float perlin_3d(float x, float y, float z) {
    int X;
    int Y;
    int Z;

    float fx = floorfrac(x, X);
    float fy = floorfrac(y, Y);
    float fz = floorfrac(z, Z);

    float u = fade(fx);
    float v = fade(fy);
    float w = fade(fz);

    using namespace lightwave::hash;

    return tri_mix(
        grad3(hash_uint3(X, Y, Z), fx, fy, fz),
        grad3(hash_uint3(X + 1, Y, Z), fx - 1.0f, fy, fz),
        grad3(hash_uint3(X, Y + 1, Z), fx, fy - 1.0f, fz),
        grad3(hash_uint3(X + 1, Y + 1, Z), fx - 1.0f, fy - 1.0f, fz),
        grad3(hash_uint3(X, Y, Z + 1), fx, fy, fz - 1.0f),
        grad3(hash_uint3(X + 1, Y, Z + 1), fx - 1.0f, fy, fz - 1.0f),
        grad3(hash_uint3(X, Y + 1, Z + 1), fx, fy - 1.0f, fz - 1.0f),
        grad3(hash_uint3(X + 1, Y + 1, Z + 1), fx - 1.0f, fy - 1.0f, fz - 1.0f),
        u,
        v,
        w);
}

/**
 * Remap the output of noise to a predictable range [-1, 1].
 * The scale values were computed experimentally by the OSL developers.
 */

inline float noise_scale1(float result) { return 0.2500f * result; }

inline float noise_scale2(float result) { return 0.6616f * result; }

inline float noise_scale3(float result) { return 0.9820f * result; }

inline float noise_scale4(float result) { return 0.8344f * result; }

/**
 * Safe Signed and Unsigned noise functions
 */

inline float snoise_1d(float p) {
    float percision_correction = 0.5f * float(fabsf(p) >= 1000000.0f);
    /**
     * Repeat Perlin noise texture every 100000.0 on each axis to prevent
     * floating point representation issues.
     *
     */
    p = fmodf(p, 100000.0f) + percision_correction;

    return noise_scale1(perlin_1d(p));
}

inline float noise_1d(float p) { return 0.5f * snoise_1d(p) + 0.5f; }

inline float snoise_2d(Point2 p) {
    Vector2 precision_correction =
        0.5f * Vector2(float(fabsf(p.x()) >= 1000000.0f),
                       float(fabsf(p.y()) >= 1000000.0f));

    /**
     * Repeat Perlin noise texture every 100000.0f on each axis to prevent
     * floating point representation issues. This causes discontinuities every
     * 100000.0f, however at such scales this usually shouldn't be noticeable.
     */
    p = (p % Point2(100000.0f)) + precision_correction;

    return noise_scale2(perlin_2d(p.x(), p.y()));
}

inline float noise_2d(Point2 p) { return 0.5f * snoise_2d(p) + 0.5f; }

inline float snoise_3d(Point p) {
    Vector precision_correction =
        0.5f * Vector(float(fabsf(p.x()) >= 1000000.0f),
                      float(fabsf(p.y()) >= 1000000.0f),
                      float(fabsf(p.z()) >= 1000000.0f));

    /**
     * Repeat Perlin noise texture every 100000.0f on each axis to prevent
     * floating point representation issues. This causes discontinuities every
     * 100000.0f, however at such scales this usually shouldn't be noticeable.
     */
    p = (p % Point(100000.0f)) + precision_correction;

    return noise_scale3(perlin_3d(p.x(), p.y(), p.z()));
}

inline float noise_3d(Point p) { return 0.5f * snoise_3d(p) + 0.5f; }

/**
 * Fractal Noise functions
 */

// TODO: fBM 1D, 3D

inline float noise_fbm(Point2 p, float detail, float roughness,
                       float lacunarity, bool normalize) {
    float fscale = 1.0f;
    float amp    = 1.0f;
    float maxamp = 0.0f;
    float sum    = 0.0f;

    for (int i = 0; i <= static_cast<int>(detail); i++) {
        float t = snoise_2d(p * fscale);
        sum += t * amp;
        maxamp += amp;
        amp *= roughness;
        fscale *= lacunarity;
    }

    float rmd = detail - floorf(detail);

    if (rmd != 0.0f) {
        float t    = snoise_2d(p * fscale);
        float sum2 = sum + t * amp;
        return normalize ? mix(0.5f * (sum / maxamp) + 0.5f,
                               0.5f * (sum2 / (maxamp + amp)) + 0.5f,
                               rmd)
                         : mix(sum, sum2, rmd);
    } else {
        return normalize ? 0.5f * (sum / maxamp) + 0.5f : sum;
    }
}

} // namespace lightwave