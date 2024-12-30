/**
 * @file transform.hpp
 * @brief Contains the Transform interface.
 */

#pragma once

#include <lightwave/color.hpp>
#include <lightwave/core.hpp>
#include <lightwave/math.hpp>

namespace lightwave {

class SimpleTransform;

/**
 * @brief Transfers points or vectors from one coordinate system to another.
 * @note This is an interface to allow time-dependent transforms (e.g., motion
 * blur), or non-linear transforms (be creative!)
 */
class Transform : public Object {
protected:
    Matrix4x4 m_transform = Matrix4x4::identity();
    Matrix4x4 m_inverse   = Matrix4x4::identity();

public:
    Transform() {}
    Transform(const Properties &) {}
    Transform(const Matrix3x3 &matrix, const bool isOrthonormal = false) {

        m_transform.setColumn(0, Vector4{ matrix.column(0), 0 });
        m_transform.setColumn(1, Vector4{ matrix.column(1), 0 });
        m_transform.setColumn(2, Vector4{ matrix.column(2), 0 });
        m_transform.setColumn(3, Vector4{ 0, 0, 0, 1 });

        Matrix3x3 inv;
        if (isOrthonormal)
            inv = matrix.transpose();
        else {
            if (auto i = invert(matrix)) {
                inv = *i;
            } else {
                lightwave_throw("transform is not invertible");
            }
        }

        m_inverse.setColumn(0, Vector4{ inv.column(0), 0 });
        m_inverse.setColumn(1, Vector4{ inv.column(1), 0 });
        m_inverse.setColumn(2, Vector4{ inv.column(2), 0 });
        m_inverse.setColumn(3, Vector4{ 0, 0, 0, 1 });
    }

    /// @brief Appends a matrix in homogeneous coordinates to this transform.
    void matrix(const Matrix4x4 &value, const bool isOrthonormal = false) {
        m_transform = value * m_transform;
        if (isOrthonormal)
            m_inverse = m_inverse * value.transpose();
        else {
            if (auto inv = invert(value)) {
                m_inverse = m_inverse * *inv;
            } else {
                lightwave_throw("transform is not invertible");
            }
        }
    }

    /// @brief Appends a translation to this transform.
    void translate(const Vector &translation) {
        // clang-format off
        m_transform = Matrix4x4 {
            1, 0, 0, translation.x(),
            0, 1, 0, translation.y(),
            0, 0, 1, translation.z(),
            0, 0, 0, 1
        } * m_transform;

        m_inverse = m_inverse * Matrix4x4 {
            1, 0, 0, -translation.x(),
            0, 1, 0, -translation.y(),
            0, 0, 1, -translation.z(),
            0, 0, 0, 1
        };
        // clang-format on
    }

    /// @brief Appends a (potentially non-uniform) scaling to this transform.
    void scale(const Vector &scaling) {
        if (scaling.product() == 0) {
            lightwave_throw("scaling is not invertible");
        }

        // clang-format off
        m_transform = Matrix4x4 {
            scaling.x(), 0, 0, 0,
            0, scaling.y(), 0, 0,
            0, 0, scaling.z(), 0,
            0, 0, 0,           1
        } * m_transform;

        m_inverse = m_inverse * Matrix4x4 {
            1 / scaling.x(), 0, 0, 0,
            0, 1 / scaling.y(), 0, 0,
            0, 0, 1 / scaling.z(), 0,
            0, 0,  0,              1
        };
        // clang-format on
    }

    /// @brief Appends a rotation around the given axis to this transform.
    void rotate(const Vector &axis, float angle) {
        const auto u    = axis.normalized();
        const float cos = std::cos(angle);
        const float sin = std::sin(angle);

        auto rotation = Matrix4x4::identity();
        for (int row = 0; row < axis.Dimension; row++) {
            for (int column = 0; column < axis.Dimension; column++) {
                rotation(row, column) =
                    (1 - cos) * u[row] * u[column] +
                    (row == column
                         ? cos
                         : (row == (column + 1) % axis.Dimension ? +1 : -1) *
                               sin *
                               u[((axis.Dimension - row) +
                                  (axis.Dimension - column)) %
                                 axis.Dimension]);
            }
        }

        m_transform = rotation * m_transform;
        m_inverse   = m_inverse * rotation.transpose();
    }

    /**
     * @brief Appends a "lookat" operation to this transform, which is useful to
     * aim cameras or light sources at other objects. The z-axis will be
     * re-oriented to be aligned with @code target - origin @endcode , and the
     * y-axis will be in the plane that the @c up vector lies in.
     */
    void lookat(const Vector &origin, const Vector &target, const Vector &up) {
        const auto direction = (target - origin).normalized();
        if (up.cross(direction).isZero()) {
            lightwave_throw(
                "lookat: direction (%s) and up vector (%s) must "
                "not be colinear",
                direction,
                up);
        }
        const auto left         = up.cross(direction).normalized();
        const auto orthogonalUp = direction.cross(left).normalized();

        Matrix4x4 matrix;
        matrix.setColumn(0, Vector4(left, 0));
        matrix.setColumn(1, Vector4(orthogonalUp, 0));
        matrix.setColumn(2, Vector4(direction, 0));
        matrix.setColumn(3, Vector4(origin, 1));

        m_transform = matrix * m_transform;

        matrix.setColumn(3, Vector4(0, 0, 0, 1));
        matrix = matrix.transpose();
        matrix.setColumn(3, Vector4(-origin, 1));

        m_inverse = m_inverse * matrix;
    }

    virtual SimpleTransform interpolate(float time) const            = 0;
    virtual const Bounds getBoundingBox(const Bounds &initial) const = 0;

    std::string toString() const override {
        return tfm::format(
            "Transform[\n"
            "  matrix = %s,\n"
            "  inverse = %s,\n"
            "]",
            indent(m_transform),
            indent(m_inverse));
    }
};

class SimpleTransform : public Object {

protected:
    Matrix4x4 m_transform;
    Matrix4x4 m_inverse;

public:
    SimpleTransform(const Matrix4x4 &transform) : m_transform(transform) {
        if (auto inv = invert(transform)) {
            m_inverse = *inv;
        } else {
            logger(EError, "transform is not invertible");
        }
    }

    /// @brief Transforms the given point.
    Point apply(const Point &point) const {
        const Vector4 result = m_transform * Vector4(Vector(point), 1);
        return Vector(result.x(), result.y(), result.z()) / result.w();
    }

    /// @brief Transforms the given vector.
    Vector apply(const Vector &vector) const {
        const Vector4 result = m_transform * Vector4(vector, 0);
        return Vector(result.x(), result.y(), result.z());
    }

    /// @brief Transforms the given normal, will not be normalized!
    Vector applyNormal(const Vector &vector) const {
        const Vector4 result = m_inverse.transpose() * Vector4(vector, 0);
        return Vector(result.x(), result.y(), result.z());
    }

    /**
     * @brief Transforms the given ray.
     * @warning The ray direction will not be normalized, as its transformed
     * length is typically useful for other tasks (e.g., instancing).
     */
    Ray apply(const Ray &ray) const {
        Ray result(ray);
        result.origin    = apply(ray.origin);
        result.direction = apply(ray.direction);
        return result;
    }

    /// @brief Applies the inverse transform to the given point.
    Point inverse(const Point &point) const {
        const Vector4 result = m_inverse * Vector4(Vector(point), 1);
        return Vector(result.x(), result.y(), result.z()) / result.w();
    }

    /// @brief Applies the inverse transform to the given vector.
    Vector inverse(const Vector &vector) const {
        const Vector4 result = m_inverse * Vector4(vector, 0);
        return Vector(result.x(), result.y(), result.z());
    }

    /**
     * @brief Transforms the given ray.
     * @warning The ray direction will not be normalized.
     */
    Ray inverse(const Ray &ray) const {
        Ray result(ray);
        result.origin    = inverse(ray.origin);
        result.direction = inverse(ray.direction);
        return result;
    }

    /// @brief Returns the determinant of this transformation.
    float determinant() const {
        return m_transform.submatrix<3, 3>(0, 0).determinant();
    }

    std::string toString() const override {
        return tfm::format(
            "SimpleTransform[\n"
            "  matrix = %s,\n"
            "  inverse = %s,\n"
            "]",
            indent(m_transform),
            indent(m_inverse));
    }
};

class StaticTransform : public Transform {

public:
    StaticTransform(const Properties &properties) : Transform(properties) {}
    StaticTransform(const Matrix3x3 &matrix, const bool isOrthonormal = false)
        : Transform(matrix, isOrthonormal) {}

    SimpleTransform interpolate(float time) const override {
        SimpleTransform T(m_transform);
        return T;
    }

    const Matrix4x4 &getTransform() const { return m_transform; }

    const Bounds getBoundingBox(const Bounds &initial) const override {

        const Bounds untransformedAABB = initial;
        if (untransformedAABB.isUnbounded()) {
            return Bounds::full();
        }

        SimpleTransform T = interpolate(0);
        Bounds result;
        for (int point = 0; point < 8; point++) {
            Point p = untransformedAABB.min();
            for (int dim = 0; dim < p.Dimension; dim++) {
                if ((point >> dim) & 1) {
                    p[dim] = untransformedAABB.max()[dim];
                }
            }
            p = T.apply(p);
            result.extend(p);
        }
        return result;
    }

    std::string toString() const override {
        return tfm::format(
            "StaticTransform[\n"
            "  matrix = %s,\n"
            "  inverse = %s,\n"
            "]",
            indent(m_transform),
            indent(m_inverse));
    }
};

class AnimatedTransform : public Transform {
protected:
    std::vector<ref<StaticTransform>> m_keyframes;

public:
    AnimatedTransform(const Properties &properties) : Transform(properties) {
        m_keyframes = properties.getChildren<StaticTransform>();
    }

    SimpleTransform interpolate(float time) const override {
        Matrix4x4 intp = m_keyframes.front()->getTransform() * (1 - time) +
                         m_keyframes.back()->getTransform() * time;
        SimpleTransform T(intp);
        return T;
    }

    const Bounds getBoundingBox(const Bounds &initial) const override {
        Bounds result;
        result.extend(m_keyframes.front()->getBoundingBox(initial));
        result.extend(m_keyframes.back()->getBoundingBox(initial));

        return result;
    }
};

} // namespace lightwave
