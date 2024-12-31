/**
 * @file math.hpp
 * @brief Contains all geometrical constructs (Points, Matrices, Rays, etc), as
 * well as commonly used mathematical constants and functions.
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>

namespace internal {

static constexpr int SnormBits = 16;

// MARK: - useful constants

/// @brief Pi
static constexpr float Pi = 3.14159265358979323846f;
/// @brief 1 / Pi
static constexpr float InvPi = 0.31830988618379067154f;
/// @brief 1 / (2 * Pi)
static constexpr float Inv2Pi = 0.15915494309189533577f;
/// @brief 1 / (4 * Pi)
static constexpr float Inv4Pi = 0.07957747154594766788f;
/// @brief Pi / 2
static constexpr float Pi2 = 1.57079632679489661923f;
/// @brief Pi / 4
static constexpr float Pi4 = 0.78539816339744830961f;
/// @brief sqrt(2)
static constexpr float Sqrt2 = 1.41421356237309504880f;

/// @brief Multiply by this constant to convert degrees into radians.
static constexpr float Deg2Rad = Pi / 180.0f;
/// @brief Multiply by this constant to convert radians into degrees.
static constexpr float Rad2Deg = 180.0f * InvPi;

static constexpr float Epsilon = 1e-4f;

/// @brief Infinity
static constexpr float Infinity = std::numeric_limits<float>::infinity();

// MARK: - utility functions

/// @brief Square root function.
inline float sqrt(float v) { return std::sqrt(v); }
/// @brief Square function, i.e., @code sqr(v) = v * v @endcode.
inline float sqr(float v) { return v * v; }

/// @brief Maximum of two integers.
inline int max(int a, int b) { return std::max(a, b); }
/// @brief Minimum of two integers.
inline int min(int a, int b) { return std::min(a, b); }

/// @brief Maximum of two floating point numbers.
inline float max(float a, float b) { return std::max(a, b); }
/// @brief Minimum of two floating point numbers.
inline float min(float a, float b) { return std::min(a, b); }

/// @brief Returns a value of @c mag with the sign of @c sgn .
inline float copysign(float mag, float sgn) { return std::copysign(mag, sgn); }
/// @brief Returns the absolute value of @c v .
inline float abs(float v) { return std::abs(v); }

/// @brief Clamps an integer @c v to lie in the interval from @c lo to @c hi .
inline int clamp(int v, int lo, int hi) { return max(lo, min(v, hi)); }
/// @brief Clamps a floating point number @c v to lie in the interval from @c lo
/// to @c hi .
inline float clamp(float v, float lo, float hi) { return max(lo, min(v, hi)); }
/// @brief Clamps a value @c v to lie in the unit interval.
inline float saturate(float v) { return clamp(v, 0.f, 1.f); }

/**
 * @brief Safe square root, which clamps negative input values to zero.
 * @note Use this when floating point errors need to be avoided, e.g., for @code
 * sin = safe_sqrt(1 - sqr(cos)) @endcode.
 */
inline float safe_sqrt(float v) { return v <= 0 ? 0 : std::sqrt(v); }
/**
 * @brief Safe arcus cosine function, which clamps the input to -1 to +1.
 * @note Use this when floating point errors need to be avoided.
 */
inline float safe_acos(float v) { return std::acos(clamp(v, -1.f, +1.f)); }

// MARK: - points and vectors

#define BUILD1(expr)                                                           \
    result;                                                                    \
    for (int i = 0; i < result.Dimension; i++)                                 \
        result[i] = expr;                                                      \
    return result;

/// @brief A point in @c D dimensions, the components of which are stored with
/// datatype @c T .
template <typename T, int D> class TPoint {
protected:
    /// @brief The components that constitute the point.
    std::array<T, D> m_data;

public:
    /// @brief The datatype used to store the components of the point.
    typedef T Type;
    /// @brief The dimensionality of the point.
    static constexpr int Dimension = D;

    /// @brief Constructs a point at the origin.
    TPoint() { std::fill(m_data.begin(), m_data.end(), Type(0)); }
    /// @brief Constructs a point from a given array.
    TPoint(const std::array<Type, Dimension> &data) : m_data(data) {}
    /// @brief Constructs a point whose components all have the value @c v .
    explicit TPoint(const Type &v) {
        std::fill(m_data.begin(), m_data.end(), v);
    }

    /// @brief Constructs a two-dimensional point.
    TPoint(const Type &x, const Type &y) : m_data({ x, y }) {}
    /// @brief Constructs a three-dimensional point.
    TPoint(const Type &x, const Type &y, const Type &z) : m_data({ x, y, z }) {}

    /// @brief Returns an array of the components of this point.
    const std::array<Type, Dimension> &data() const { return m_data; }
    /// @brief Returns an array of the components of this point that can be
    /// modified.
    std::array<Type, Dimension> &data() { return m_data; }

    /// @brief Access a component of this point, with an index ranging from @c 0
    /// to @code Dimension - 1 @endcode .
    const Type &operator[](int i) const { return m_data[i]; }
    /// @brief Access a component of this point that can be modified, with an
    /// index ranging from @c 0 to @code Dimension - 1 @endcode .
    Type &operator[](int i) { return m_data[i]; }

    const Type &x() const {
        static_assert(Dimension >= 1);
        return m_data[0];
    }
    const Type &y() const {
        static_assert(Dimension >= 2);
        return m_data[1];
    }
    const Type &z() const {
        static_assert(Dimension >= 3);
        return m_data[2];
    }
    const Type &w() const {
        static_assert(Dimension >= 4);
        return m_data[3];
    }

    Type &x() {
        static_assert(Dimension >= 1);
        return m_data[0];
    }
    Type &y() {
        static_assert(Dimension >= 2);
        return m_data[1];
    }
    Type &z() {
        static_assert(Dimension >= 3);
        return m_data[2];
    }
    Type &w() {
        static_assert(Dimension >= 4);
        return m_data[3];
    }

    /// @brief Returns the elementwise minimum of two points.
    friend auto elementwiseMin(const TPoint &a, const TPoint &b) {
        TPoint BUILD1(std::min(a[i], b[i]))
    }
    /// @brief Returns the elementwise maximum of two points.
    friend auto elementwiseMax(const TPoint &a, const TPoint &b) {
        TPoint BUILD1(std::max(a[i], b[i]))
    }

    /// @brief Checks whether two points are exactly identical.
    bool operator==(const TPoint &other) const {
        return m_data == other.m_data;
    }
    /// @brief Checks whether two points are not exactly identical.
    bool operator!=(const TPoint &other) const {
        return m_data != other.m_data;
    }

    /// @brief Returns whether the point lies at the origin, i.e., all
    /// components are zero.
    bool isZero() const {
        return std::all_of(
            m_data.begin(), m_data.end(), [](const Type &v) { return v == 0; });
    }
};

/// @brief A vector in @c D dimensions, the components of which are stored with
/// datatype @c T .
template <typename Type, int Dimension>
class TVector : public TPoint<Type, Dimension> {
protected:
    using TPoint<Type, Dimension>::m_data;

public:
    using TPoint<Type, Dimension>::TPoint;

    TVector() : TPoint<Type, Dimension>(0) {}

    explicit TVector(const TPoint<Type, Dimension> &point)
        : TPoint<Type, Dimension>(point.data()) {}

    /// @brief Computes the dot product (aka scalar product) with respect to
    /// another vector.
    Type dot(const TVector &other) const {
        Type result = 0;
        for (int i = 0; i < Dimension; i++)
            result += m_data[i] * other.m_data[i];
        return result;
    }

    /// @brief Computes the cross product with respect to another vector.
    TVector cross(const TVector &other) const {
        static_assert(Dimension == 3);
        return { this->y() * other.z() - this->z() * other.y(),
                 this->z() * other.x() - this->x() * other.z(),
                 this->x() * other.y() - this->y() * other.x() };
    }

    /// @brief Computes the squared length of this vector.
    Type lengthSquared() const { return dot(*this); }
    /// @brief Computes the length of this vector.
    Type length() const { return sqrt(lengthSquared()); }
    /// @brief Returns a normalized copy of this vector.
    auto normalized() const { return *this / length(); }

    friend auto operator-(const TVector &a) { TVector BUILD1(-a[i]) }
    friend auto operator*(const Type &a, const TVector &b) {
        TVector BUILD1(a * b[i])
    }
    friend auto operator*(const TVector &a, Type b) { TVector BUILD1(a[i] * b) }
    friend auto operator/(const TVector &a, Type b) { TVector BUILD1(a[i] / b) }
    friend auto operator+(const TVector &a, const TVector &b) {
        TVector BUILD1(a[i] + b[i])
    }
    friend auto operator-(const TVector &a, const TVector &b) {
        TVector BUILD1(a[i] - b[i])
    }
    friend auto operator*(const TVector &a, const TVector &b) {
        TVector BUILD1(a[i] * b[i])
    }
    friend auto operator/(const TVector &a, const TVector &b) {
        TVector BUILD1(a[i] / b[i])
    }

    auto operator*=(const Type &other) { return *this = *this * other; }
    auto operator/=(const Type &other) { return *this = *this / other; }
    auto operator+=(const TVector &other) { return *this = *this + other; }
    auto operator-=(const TVector &other) { return *this = *this - other; }
    auto operator*=(const TVector &other) { return *this = *this * other; }
    auto operator/=(const TVector &other) { return *this = *this / other; }

    /// @brief Converts the vector to a different datatype.
    template <typename Type2> auto cast() const {
        TVector<Type2, Dimension> BUILD1(Type2(m_data[i]))
    }
};

#undef BUILD_VECTOR

template <typename Type, int Dimension>
auto operator+(const TPoint<Type, Dimension> &a,
               const TVector<Type, Dimension> &b) {
    TPoint<Type, Dimension> BUILD1(a[i] + b[i])
}

template <typename Type, int Dimension>
auto operator+=(TPoint<Type, Dimension> &a, const TVector<Type, Dimension> &b) {
    a = a + b;
}

template <typename Type, int Dimension>
auto operator-(const TPoint<Type, Dimension> &a,
               const TVector<Type, Dimension> &b) {
    TPoint<Type, Dimension> BUILD1(a[i] - b[i])
}

template <typename Type, int Dimension>
auto operator-=(TPoint<Type, Dimension> &a, const TVector<Type, Dimension> &b) {
    a = a - b;
}

template <typename Type, int Dimension>
auto operator-(const TPoint<Type, Dimension> &a,
               const TPoint<Type, Dimension> &b) {
    TVector<Type, Dimension> BUILD1(a[i] - b[i])
}

struct Snorm {
    int16_t m_bits;

    Snorm() : m_bits(0) {}

    Snorm(int bits) : m_bits(static_cast<int16_t>(bits)) {}
    Snorm(float value) {
        m_bits = (int16_t) round(value * ((1 << (SnormBits - 1)) - 1));
    }

    int16_t data() const { return m_bits; }

    static Snorm flooredSnorm(float value) {
        return Snorm((int16_t) floor(value * ((1 << (SnormBits - 1)) - 1)));
    }

    operator float() const {
        return clamp((float) m_bits / ((1 << (SnormBits - 1)) - 1), -1.0, 1.0);
    }
};

/// @brief A three-dimensional point with floating point components.
using Point = TPoint<float, 3>;

/// @brief A two-dimensional vector with floating point components.
using Vector2 = TVector<float, 2>;
/// @brief A two-dimensional vector with floating point components quantized as
/// 16-bit integer.
using Vector2s = TVector<Snorm, 2>;
/// @brief A three-dimensional vector with floating point components.
using Vector = TVector<float, 3>;
/// @brief A three-dimensional vector with integer components.
using Vector3i = TVector<int, 3>;

struct Bounds {
    Point m_min;
    Point m_max;

    Bounds() : m_min(0), m_max(0) {};
    Bounds(Point min, Point max) : m_min(min), m_max(max) {};

    /// @brief Extends this bounding box to also cover the region of another
    /// bounding box.
    void extend(const Point &p) {
        m_min = elementwiseMin(m_min, p);
        m_max = elementwiseMax(m_max, p);
    }

    /// @brief Returns the diagonal of the bounding box, i.e., @code max - min
    /// @endcode .
    Vector diagonal() const { return m_max - m_min; }

    /// @brief Returns the lower corner of this bounding box.
    const Point &min() const { return m_min; }
    /// @brief Returns the upper corner of this bounding box.
    const Point &max() const { return m_max; }
    /// @brief Returns a reference to the lower corner of this bounding box.
    Point &min() { return m_min; }
    /// @brief Returns a reference to the upper corner of this bounding box.
    Point &max() { return m_max; }
};

} // namespace internal
