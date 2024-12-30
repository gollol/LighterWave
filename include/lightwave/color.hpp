/**
 * @file color.hpp
 * @brief Contains the Color class, used to represent and work with RGB colors.
 */

#pragma once

#include <lightwave/math.hpp>

#include <ostream>
#include <array>

namespace lightwave {

#define BUILD(expr)                                                            \
    Color result;                                                              \
    for (int i = 0; i < result.NumComponents; i++)                             \
        result[i] = expr;                                                      \
    return result;

/// @brief Represents RGB colors in linear color space.
class Color {
public:
    static constexpr int NumComponents = 3;

private:
    /// @brief Contains the RGB components of this color, in that sequence.
    std::array<float, NumComponents> m_data;

public:
    /// @brief Create black color.
    Color() { std::fill(m_data.begin(), m_data.end(), 0.0f); }
    /// @brief Create gray color with brightness @c v .
    explicit Color(float v) { std::fill(m_data.begin(), m_data.end(), v); }

    /// @brief Create color with the provided RGB values.
    Color(float r, float g, float b) : m_data({ r, g, b }) {}
    /// @brief Interpret a @ref Vector as color ( @c x corresponds to @c r , @c
    /// y to @c g , and @c z to @c b ).
    explicit Color(const Vector &vec) : m_data(vec.data()) {}

    /// @brief Returns an array of the RGB values of this color, in that
    /// sequence.
    const std::array<float, 3> &data() const { return m_data; }
    /// @brief Returns an array of the RGB values of this color that can be
    /// modified, in that sequence.
    std::array<float, 3> &data() { return m_data; }

    /// @brief Access a component of this color, with an index either 0 (red), 1
    /// (green), or 2 (blue).
    const float &operator[](int i) const { return m_data[i]; }
    /// @brief Access a component of this color that can be modified, with an
    /// index either 0 (red), 1 (green), or 2 (blue).
    float &operator[](int i) { return m_data[i]; }

    /// @brief Get the red component of this color.
    const float &r() const { return m_data[0]; }
    /// @brief Get the green component of this color.
    const float &g() const { return m_data[1]; }
    /// @brief Get the blue component of this color.
    const float &b() const { return m_data[2]; }

    /// @brief Get the red component of this color that can be modified.
    float &r() { return m_data[0]; }
    /// @brief Get the green component of this color that can be modified.
    float &g() { return m_data[1]; }
    /// @brief Get the blue component of this color that can be modified.
    float &b() { return m_data[2]; }

    /// @brief Multiply all color components by a given scalar.
    friend Color operator*(float a, const Color &b) { BUILD(a * b[i]) }
    /// @brief Multiply all color components by a given scalar.
    friend Color operator*(const Color &a, float b) { BUILD(a[i] * b) }
    /// @brief Divide all color components by a given scalar.
    friend Color operator/(const Color &a, float b) { BUILD(a[i] / b) }
    /// @brief Add two colors component-wise.
    friend Color operator+(const Color &a, const Color &b) { BUILD(a[i] + b[i]) }
    /// @brief Subtract two colors component-wise.
    friend Color operator-(const Color &a, const Color &b) { BUILD(a[i] - b[i]) }
    /// @brief Multiply two colors component-wise.
    friend Color operator*(const Color &a, const Color &b) { BUILD(a[i] * b[i]) }
    /// @brief Divide two colors component-wise.
    friend Color operator/(const Color &a, const Color &b) { BUILD(a[i] / b[i]) }

    /// @brief Returns the component-wise maximum of two colors.
    friend Color max(const Color &a, const Color &b) {
        BUILD(std::max(a[i], b[i]))
    }
    /// @brief Returns the component-wise minimum of two colors.
    friend Color min(const Color &a, const Color &b) {
        BUILD(std::min(a[i], b[i]))
    }

    /// @brief Clamps each component to lie in the range 0 to 1.
    friend Color saturate(const Color &a) { BUILD(saturate(a[i])) }

    /// @brief Multiplies the color by a scalar.
    auto operator*=(const float &other) { return *this = *this * other; }
    /// @brief Divides the color by a scalar.
    auto operator/=(const float &other) { return *this = *this / other; }
    /// @brief Adds another color component-wise.
    auto operator+=(const Color &other) { return *this = *this + other; }
    /// @brief Subtracts another color component-wise.
    auto operator-=(const Color &other) { return *this = *this - other; }
    /// @brief Multiplies another color component-wise.
    auto operator*=(const Color &other) { return *this = *this * other; }
    /// @brief Divides another color component-wise.
    auto operator/=(const Color &other) { return *this = *this / other; }

    /// @brief Checks whether two colors are exactly identical.
    bool operator==(const Color &other) const { return m_data == other.m_data; }
    /// @brief Checks whether two colors are not exactly identical.
    bool operator!=(const Color &other) const { return m_data != other.m_data; }

    /// @brief Returns the luminance of this color.
    float luminance() const {
        return r() * 0.212671f + g() * 0.715160f + b() * 0.072169f;
    }

    bool isValid() const {
        for (int i = 0; i < NumComponents; i++)
            if (!std::isfinite(m_data[i]) || m_data[i] < 0)
                return false;
        return true;
    }

    /// @brief Returns the arithmetic mean of the components of this color.
    float mean() const { return (1 / 3.f) * (r() + g() + b()); }

    /// @brief Creates black color (i.e., all components 0).
    static Color black() { return Color(0); }
    /// @brief Creates white color (i.e., all components 1).
    static Color white() { return Color(1); }
};

/// @brief Linearly blend between two colors @c a (t = 0) and @c b (t = 1).
static Color lerp(const Color &a, const Color &b, float t) {
    return a * (1 - t) + b * t;
}

#undef BUILD

/// @brief Print a given color to an output stream.
static std::ostream &operator<<(std::ostream &os, const Color &color) {
    os << "Color[" << color.r() << ", " << color.g() << ", " << color.b()
       << "]";
    return os;
}

/// @brief Converrts a given color value in RGB to its corresponding HSV representation
/// @param color The original color in the RGB space
/// @return The converted color in the HSV space
/// @note See https://de.wikipedia.org/wiki/HSV-Farbraum#Umrechnung_RGB_in_HSV/HSL
/// @note Also see https://github.com/iRath96/raymond/blob/fd6341d2c1c63698705343c3aabfffb1c13905a5/raymond/device/utils/color.hpp#L9
inline Color convertRGBToHSV(const Color &color) {
    const float MAX = max(color.r(), max(color.g(), color.b()));
    const float MIN = min(color.r(), min(color.g(), color.b()));
    const float DELTA = MAX - MIN;

    float H;
    float S;
    float V = MAX;

    // Get value for H
    if (MAX == MIN) {
        H = 0.0f;
    } else if (MAX == color.r()) {
        H = 0 + ((color.g() - color.b()) / DELTA);
    } else if (MAX == color.g()) {
        H = 2 + ((color.b() - color.r()) / DELTA);
    } else if (MAX == color.b()) {
        H = 4 + ((color.r() - color.g()) / DELTA);
    } else {
        H = 0.0f;
        logger(EWarn, "[convertRGBToHSV] MAX (%s) is not either of the RGB values! How did we get here?", color);
    }

    // Scale H value to be between 0.0f and 1.0f
    H = H / 6.0f;

    // Get value for S
    if (MAX == MIN) {
        S = 0.0f;
    } else {
        S = DELTA / MAX;
    }

    return Color(H, S, V);
}

/// @brief Converrts a given color value in HSV to its corresponding RGB representation
/// @param color The original color in the HSV space
/// @return The converted color in the RGB space
/// @note See https://de.wikipedia.org/wiki/HSV-Farbraum#Umrechnung_RGB_in_HSV/HSL
/// @note Also see https://github.com/iRath96/raymond/blob/fd6341d2c1c63698705343c3aabfffb1c13905a5/raymond/device/utils/color.hpp#L43
inline Color convertHSVToRGB(const Color &color) {
    // To start, H is in the range of [0,1]
    const float H = color.r() * 6.0f;
    const float S = color.g();
    const float V = color.b();

    if (S == 0.0f) {
        return Color(V, V, V);
    }

    const int h_i = static_cast<int>(floorf(H));
    const float f = H - h_i;

    const float p = V * (1 - S);
    const float q = V * (1 - S * f);
    const float t = V * (1 - S * (1 - f));

    switch (h_i) {
    case 0:
    case 6:
        return Color(V, t, p);
    case 1:
        return Color(q, V, p);
    case 2:
        return Color(p, V, t);
    case 3:
        return Color(p, q, V);
    case 4:
        return Color(t, p, V);
    case 5:
        return Color(V, p, q);
    
    default:
        logger(EWarn, "[convertRGBToHSV] Can't convert HSV color to RGB: %s", color);
        return Color::black();
    }
}

} // namespace lightwave

namespace std {

/// @brief Checks whether all components of the color are finite.
static bool isfinite(const lightwave::Color &c) {
    for (int chan = 0; chan < c.NumComponents; chan++)
        if (!std::isfinite(c[chan]))
            return false;
    return true;
}

} // namespace std
