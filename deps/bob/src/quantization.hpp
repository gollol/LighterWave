#pragma once

#include "math.hpp"
#include <fstream>

namespace internal {

// TODO: do this automatically, maybe?
inline Vector2s encode(const Vector2 &w) {
    return { Snorm(w.x()), Snorm(w.y()) };
}

class Spherical {
    Vector2s m_encoded;

public:
    Spherical(const Vector &v) {
        m_encoded = encode(Vector2((acos(v.z()) * InvPi - 0.5f) * 2.0f,
                                   atan2(v.y(), v.x()) * InvPi));
    }

    operator Vector() const {
        float theta = (float(m_encoded.x()) * 0.5f + 0.5f) * Pi;
        float phi   = float(m_encoded.y()) * Pi;

        return { static_cast<float>(sin(theta) * cos(phi)),
                 static_cast<float>(sin(theta) * sin(phi)),
                 static_cast<float>(cos(theta)) };
    }
};

class Oct32 {
protected:
    Vector2s m_encoded;

public:
    Oct32() = default;
    Oct32(const Vector &v) {
        float c = 1.0f / (abs(v.x()) + abs(v.y()) + abs(v.z()));

        if (v.z() >= 0.0f) {
            m_encoded = encode(Vector2(v.x() * c, v.y() * c));
        } else {
            m_encoded = encode(Vector2(copysign(1.0f - abs(v.y()) * c, v.x()),
                                       copysign(1.0f - abs(v.x()) * c, v.y())));
        }
    }

    Oct32(const Vector2s &v) : m_encoded(v) {}

    operator Vector() const {
        Vector result;
        result.x() = float(m_encoded.x());
        result.y() = float(m_encoded.y());
        result.z() = 1.0f - (abs(result.x()) + abs(result.y()));

        if (result.z() < 0.0f) {
            float oldX = result.x();
            result.x() = copysign(1.0f - abs(result.y()), oldX);
            result.y() = copysign(1.0f - abs(oldX), result.y());
        }
        return result.normalized();
    }
};

class Oct32p {
    Vector2s m_encoded;

public:
    Oct32p(const Vector &v) {
        Vector2s projected;
        float c = 1.0f / (abs(v.x()) + abs(v.y()) + abs(v.z()));

        if (v.z() >= 0.0f) {
            projected = { Snorm::flooredSnorm(v.x() * c),
                          Snorm::flooredSnorm(v.y() * c) };
        } else {
            projected = {
                Snorm::flooredSnorm(copysign(1.0f - abs(v.x()) * c, v.x())),
                Snorm::flooredSnorm(copysign(1.0f - abs(v.y()) * c, v.y()))
            };
        }

        m_encoded      = projected;
        float accuracy = 0.0f;

        for (int16_t i = 0; i < 2; i++) {
            for (int16_t j = 0; j < 2; j++) {
                Vector2s current = projected;
                current.x().m_bits += i;
                current.y().m_bits += j;

                Vector decoded         = Oct32(current);
                float current_accuracy = abs(decoded.dot(v));
                if (current_accuracy > accuracy) {
                    accuracy  = current_accuracy;
                    m_encoded = current;
                }
            }
        }
    }

    operator Vector() const {
        Vector result;
        result.x() = float(m_encoded.x());
        result.y() = float(m_encoded.y());
        result.z() = 1.0f - (abs(result.x()) + abs(result.y()));

        if (result.z() < 0.0f) {
            result.x() = copysign(1.0f - abs(result.x()), result.x());
            result.y() = copysign(1.0f - abs(result.y()), result.y());
        }
        return result.normalized();
    }
};

struct Vertex {
    Point position;
    Vector2 uv;
    Oct32 normal;

    Vertex() = default;

    Vertex(const Point &p, const Vector2 &uv, const Vector &n)
        : position(p), uv(uv), normal(Oct32(n)) {}

    void write(std::ofstream &file) {
        file.write((char *) &position, sizeof(position));
        file.write((char *) &uv, sizeof(uv));
        file.write((char *) &normal, sizeof(normal));
    }
};

} // namespace internal