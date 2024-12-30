/**
 * @file hash.hpp
 * @brief Includes hashing functions used by data structures or to seed random
 * number generators.
 */

#pragma once

#include <lightwave/core.hpp>

namespace lightwave::hash {

/// @brief The FNV-1a (Fowler–Noll–Vo) hash function with 64-bit state.
class fnv1a {
    /// @brief The current state of the hash function.
    uint64_t hash = 0xCBF29CE484222325;

public:
    /// @brief Constructs a FNV-1a hash for a given sequence of integers.
    template <typename... Ts> fnv1a(Ts... args) { (operator<<(args), ...); }
    
    /// @brief Updates the state by hashing a given integer.
    template <typename T> fnv1a &operator<<(T data) {
        for (size_t pos = 0; pos < sizeof(T); pos++) {
            uint8_t byte = data >> (8 * pos);
            hash         = (hash ^ byte) * 0x100000001b3;
        }
        return *this;
    }

    /// @brief Returns the current state of the hash function.
    operator uint64_t() { return hash; }
};


/// @note The following code is adapted from https://projects.blender.org/blender/blender/src/commit/4439e395301672e7cb71ca15974c0d9aef294e6e/intern/cycles/util/hash.h

/**
 * Jenkins Lookup3 Hash Functions
 * Source: http://burtleburtle.net/bob/c/lookup3.c
*/
constexpr static const float fac16incl = 1.f/(float)UINT16_MAX;
constexpr static const float fac16excl = 1.f/65536.0f;

constexpr static const float fac32incl = 1.0f/(float)0xFFFFFFFFu;
constexpr static const float fac32excl = 1.0f/4294967808.0f;

#define rot(x, k) (((x) << (k)) | ((x) >> (32 - (k))))

#define mix(a, b, c) \
  { \
    a -= c; \
    a ^= rot(c, 4); \
    c += b; \
    b -= a; \
    b ^= rot(a, 6); \
    a += c; \
    c -= b; \
    c ^= rot(b, 8); \
    b += a; \
    a -= c; \
    a ^= rot(c, 16); \
    c += b; \
    b -= a; \
    b ^= rot(a, 19); \
    a += c; \
    c -= b; \
    c ^= rot(b, 4); \
    b += a; \
  } \
  ((void)0)

#define hashfinal(a, b, c) \
  { \
    c ^= b; \
    c -= rot(b, 14); \
    a ^= c; \
    a -= rot(c, 11); \
    b ^= a; \
    b -= rot(a, 25); \
    c ^= b; \
    c -= rot(b, 16); \
    a ^= c; \
    a -= rot(c, 4); \
    b ^= a; \
    b -= rot(a, 14); \
    c ^= b; \
    c -= rot(b, 24); \
  } \
  ((void)0)

inline unsigned int hash_uint(unsigned int kx) {
    unsigned int a, b, c;
    a = b = c = 0xdeadbeef + (1 << 2) + 13;

    a += kx;
    hashfinal(a, b, c);

    return c;
}

inline unsigned int hash_uint2(unsigned int kx, unsigned int ky) {
    unsigned int a, b, c;
    a = b = c = 0xdeadbeef + (2 << 2) + 13;

    b += ky;
    a += kx;
    hashfinal(a, b, c);

    return c;
}

inline unsigned int hash_uint3(unsigned int kx, unsigned int ky, unsigned int kz) {
    unsigned int a, b, c;
    a = b = c = 0xdeadbeef + (3 << 2) + 13;

    c += kz;
    b += ky;
    a += kx;
    hashfinal(a, b, c);

    return c;
}

#undef rot
#undef hashfinal
#undef mix


/* [0, uint_max] -> [0.0, 1.0) */
inline float uint_to_float_excl(unsigned int n) {
    // Note: we divide by 4294967808 instead of 2^32 because the latter
    // leads to a [0.0, 1.0] mapping instead of [0.0, 1.0) due to floating
    // point rounding error. 4294967808 unfortunately leaves (precisely)
    // one unused ulp between the max number this outputs and 1.0, but
    // that's the best you can do with this construction.
    return (float)n * fac32excl;
}

/* [0, uint_max] -> [0.0, 1.0] */
inline float uint_to_float_incl(unsigned int n) {
    return (float)n * fac32incl;
}

/**
 * Utility conversion functions from Blender
*/

inline int __float_as_int(float f) {
    union {
        int i;
        float f;
    } u;
    u.f = f;
    return u.i;
}

inline float __int_as_float(int i) {
    union {
        int i;
        float f;
    } u;
    u.i = i;
    return u.f;
}

inline unsigned int __float_as_uint(float f) {
    union {
        unsigned int i;
        float f;
    } u;
    u.f = f;
    return u.i;
}

inline float __uint_as_float(unsigned int i) {
    union {
        unsigned int i;
        float f;
    } u;
    u.i = i;
    return u.f;
}

/**
 * Hashing unsigned int or unsigned int[234] into a float in the range [0, 1].
*/

inline float hash_uint_to_float(unsigned int kx) {
  return uint_to_float_incl(hash_uint(kx));
}

inline float hash_uint2_to_float(unsigned int kx, unsigned int ky) {
  return uint_to_float_incl(hash_uint2(kx, ky));
}

inline float hash_uint3_to_float(unsigned int kx, unsigned int ky, unsigned int kz) {
  return uint_to_float_incl(hash_uint3(kx, ky, kz));
}

/**
 * Hashing float or float[234] into a float in the range [0, 1].
*/

inline float hash_float_to_float(float k) {
  return hash_uint_to_float(__float_as_uint(k));
}

inline float hash_Point2_to_float(Point2 k) {
  return hash_uint2_to_float(__float_as_uint(k.x()), __float_as_uint(k.y()));
}

inline float hash_Point_to_float(Point k) {
  return hash_uint3_to_float(__float_as_uint(k.x()), __float_as_uint(k.y()), __float_as_uint(k.z()));
}


} // namespace lightwave::hash
