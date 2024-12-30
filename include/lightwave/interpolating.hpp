/**
 * @file interpolating.hpp
 * @brief Contains some interpolation functions
 */

#pragma once

#include <lightwave/math.hpp>

namespace lightwave {

    #define BUILD1(expr)                                                       \
    result;                                                                    \
    for (int i = 0; i < result.Dimension; i++)                                 \
        result[i] = expr;                                                      \
    return result;

    template<typename T, typename Type, int D>
    TPoint<T, D> operator*(const Type &a, const TPoint<T, D> &b) {TPoint<T, D> BUILD1(a * b[i])}

    template<typename T, typename Type, int D>
    TPoint<T, D> operator*(const TPoint<T, D> &a, Type b) {TPoint<T, D> BUILD1(a[i] * b)}

    template<typename T, int D>
    TPoint<T, D> operator+(const TPoint<T, D> &a, const TPoint<T, D> &b) {TPoint<T, D> BUILD1(a[i] + b[i])}

    template<typename T, int D>
    TPoint<T, D> operator-(const TPoint<T, D> &a) {TPoint<T, D> BUILD1(-a[i])}

    #undef BUILD1
    

    template<typename T>
    inline T interpolateCubicBSpline(const float t, const std::array<T, 4> &values) {
        const T b0 = values[0];
        const T b1 = values[1];
        const T b2 = values[2];
        const T b3 = values[3];
        const float t2 = sqr(t);
        const float t3 = t*t2;
        return 1.f/6.f * (-t3*(b0-3.f*b1+3.f*b2-b3) + 3.f*t2*(b0-2.f*b1+b2) - 3.f*t*(b0-b2) + (b0 +4.f*b1+b2));
    }

    template<typename T>
    inline T interpolateBezierCurve(const float t, const std::array<T, 4> &values) {
        const T b0 = values[0];
        const T b1 = values[1];
        const T b2 = values[2];
        const T b3 = values[3];
        const float t2 = sqr(t);
        return b0 - 3.f*t*(b0-b1) + 3.f*t2*(b0-2.f*b1+b2) - t*t2*(b0-3.f*b1+3.f*b2-b3);
    }

    template<typename T>
    inline T tangentToBezierCurve(const float t, const std::array<T, 4> &values) {
        const T b0 = values[0];
        const T b1 = values[1];
        const T b2 = values[2];
        const T b3 = values[3];
        return - 3.f*(b0-b1) + 6.f*t*(b0-2.f*b1+b2) - 3.f*t*t*(b0-3.f*b1+3.f*b2-b3);
    }

    template<typename T>
    inline T normalToBezierCurve(const float t, const std::array<T, 4> &values) {
        const T b0 = values[0];
        const T b1 = values[1];
        const T b2 = values[2];
        const T b3 = values[3];
        return 6.f*(b0-2.f*b1+b2) - 6.f*t*(b0-3.f*b1+3.f*b2-b3);
    }

    /// @brief Cardinal interpolation where we require values[3] to be the "last"/"rightmost" element and values[0] to be the "first"/"leftmost" one.
    /// If positions is one-dimensional, use Floats4!
    /// @tparam T The type of the to be interpolated values (e.g. Color, Point)
    /// @param t How far we are between the control points
    /// @param values The values to be interpolated (also called control points)
    /// @param positions The coordinates of the control points / values
    /// @param c Determines the shape of the spline, usually taken to be a half.
    /// @return 
    template<typename T, int D>
    inline T interpolateCardinal(const float t, const std::array<T, 4> &values, const std::array<TPoint<float, D>, 4> &positions, const float c) {
        if(positions[1] == positions[2])
            return 0.5f*(values[1]+values[2]);
        if(positions[0] == positions[2] || positions[1] == positions[3])
            lightwave_throw("Cardinal interpolation requires pos[0] != pos[2] and pos[1] != pos[3]");
        const float h00 = (1.f + 2.f*t)*sqr(1.f-t);
        const float h10 = t*sqr(1.f-t);
        const float h01 = sqr(t)*(3.f-2.f*t);
        const float h11 = sqr(t)*(t-1.f);

        const float cInv = 1-c;
        const float dist = (positions[2]-positions[1]).length();

        const auto m1 = cInv*(values[2]-values[0])/(positions[2]-positions[0]).length();
        const auto m2 = cInv*(values[3]-values[1])/(positions[3]-positions[1]).length();

        return h00 * values[1] + h10 * dist * m1 + h01 * values[2] + h11 * dist * m2;
    }

    template<int D>
    inline TPoint<float, D> interpolateCardinalCurve(const float t, const std::array<TPoint<float, D>, 4> &controlPoints, const float c) {
        return interpolateCardinal(t, controlPoints, controlPoints, c);
    }

    template<typename T, int D>
    inline T interpolateCatmullRom(const float t, const std::array<T, 4> &values, const std::array<TPoint<float, D>, 4> &positions) {
        return interpolateCardinal(t, values, positions, 0.5f);
    }

    template<int D>
    inline TPoint<float, D> interpolateCatmullRomCurve(const float t, const std::array<TPoint<float, D>, 4> &controlPoints) {
        return interpolateCardinal(t, controlPoints, controlPoints, 0.5f);
    }

    template<typename T>
    inline T interpolateLinear(float t, const T start, const T end) {
        return start + t * (end - start);
    }
    template<typename T>
    inline T interpolateLinearPrime(float t, const T start, const T end) {
        return end - start;
    }

    template<typename T>
    inline T interpolateParabolic(float t, const T start, const T end) {
        t = 4 * (1 - t) * t;
        return start + t * (end - start);
    }
    template<typename T>
    inline T interpolateParabolicPrime(float t, const T start, const T end) {
        t = 4 - 8 * t;
        return t * (end - start);
    }

    template<typename T>
    inline T interpolateCubic(float t, const T start, const T end) {
        float t2 = t * t;
        t = 0.038 + 8.615 * t - 23.14 * t2 + 15.46 * t * t2;
        return start + t * (end - start);
    }
    template<typename T>
    inline T interpolateCubicPrime(float t, const T start, const T end) {
        t = 8.615 - 46.28 * t + 46.38 * t * t;
        return t * (end - start);
    }

    template<typename T>
    inline T interpolateEase(float t, const T start, const T end) {
        t = t * t * (3.f - 2.f * t);
        return start + t * (end - start);
    }
    template<typename T>
    inline T interpolateEasePrime(float t, const T start, const T end) {
        t = t * (6.f - 6.f * t);
        return t * (end - start);
    }

    template<typename T>
    inline T interpolateWave(float t, const int halfpeaks, const T start, const T end) {
        t = 0.5 - 0.5 * cosf(halfpeaks * Pi * t);
        return start + t * (end - start);
    }
    template<typename T>
    inline T interpolateWavePrime(float t, const int halfpeaks, const T start, const T end) {
        t = 0.5 * sinf(halfpeaks * Pi * t) * halfpeaks * Pi;
        return t * (end - start);
    }

}  // namespace lightwave