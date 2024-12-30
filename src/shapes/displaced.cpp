#include <lightwave.hpp>

#include "../core/plyparser.hpp"
#include "accel.hpp"
#include "bob.hpp"

#include <float.h>

//-------------------------------------------------------------------------
// This file contains the implementation for Tessellation-free displacement
// mapping
// Original Paper: https://perso.telecom-paristech.fr/boubek/papers/TFDM/
// Different implementation: https://github.com/shocker-0x15/GfxExp/tree/master
//-------------------------------------------------------------------------

namespace lightwave {

using Matrix2x2 = TMatrix<float, 2, 2>;

// MARK: - Utility functions

/// @brief Returns vector with entry wise abs values of input vector
Vector abs(const Vector &vec) {
    return Vector(abs(vec[0]), abs(vec[1]), abs(vec[2]));
}

float Sign(lightwave::Vector2 p1, lightwave::Vector2 p2,
           lightwave::Vector2 p3) {
    return (p1.x() - p3.x()) * (p2.y() - p3.y()) -
           (p2.x() - p3.x()) * (p1.y() - p3.y());
}

enum CollisionType {
    SquareOutsideTriangle,
    SquareInsideTriangle,
    SquareOverlappingTriangle
};

bool PointInTriangle(lightwave::Vector2 pt, lightwave::Vector2 v1,
                     lightwave::Vector2 v2, lightwave::Vector2 v3) {
    float d1, d2, d3;
    bool has_neg, has_pos;

    d1 = Sign(pt, v1, v2);
    d2 = Sign(pt, v2, v3);
    d3 = Sign(pt, v3, v1);

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}

bool Contains(lightwave::Vector2 p1, lightwave::Vector2 p2,
              lightwave::Vector2 p3, float w) {
    return PointInTriangle(lightwave::Vector2(-w, -w), p1, p2, p3) &&
           PointInTriangle(lightwave::Vector2(w, -w), p1, p2, p3) &&
           PointInTriangle(lightwave::Vector2(-w, w), p1, p2, p3) &&
           PointInTriangle(lightwave::Vector2(w, w), p1, p2, p3);
}

/// @brief Returns if and how a quad defined by its center point c and width w
/// collides with the 2D triangle defined by the positions p1, p2, p3 and the
/// normals n1, n2, n3
CollisionType Collision(lightwave::Vector2 p1, lightwave::Vector2 p2,
                        lightwave::Vector2 p3, lightwave::Vector2 n1,
                        lightwave::Vector2 n2, lightwave::Vector2 n3,
                        lightwave::Vector2 c, float w) {
    p1 -= c;
    p2 -= c;
    p3 -= c;

    bool b1 = std::min(w, std::max(p3.x(), std::max(p1.x(), p2.x()))) <=
              std::max(-w, std::min(p1.x(), std::min(p2.x(), p3.x())));
    bool b2 = std::min(w, std::max(p3.y(), std::max(p1.y(), p2.y()))) <=
              std::max(-w, std::min(p1.y(), std::min(p2.y(), p3.y())));

    float d1 = n1.dot(p1 + lightwave::Vector2(n1.x() >= 0.0f ? w : -w,
                                              n1.y() >= 0.0f ? w : -w));
    float d2 = n2.dot(p2 + lightwave::Vector2(n2.x() >= 0.0f ? w : -w,
                                              n2.y() >= 0.0f ? w : -w));
    float d3 = n3.dot(p3 + lightwave::Vector2(n3.x() >= 0.0f ? w : -w,
                                              n3.y() >= 0.0f ? w : -w));

    if (b1 || b2 || d1 <= 0.0f || d2 <= 0.0f || d3 <= 0.0f)
        return SquareOutsideTriangle;

    if (Contains(p1, p2, p3, w))
        return SquareInsideTriangle;

    return SquareOverlappingTriangle;
}

// MARK: - Affine Arithmetic stuff

struct AffineVector {
    lightwave::Vector x_c;
    lightwave::Vector x_u;
    lightwave::Vector x_v;
    lightwave::Vector x_w;

    AffineVector() {
        x_c = lightwave::Vector();
        x_u = lightwave::Vector();
        x_v = lightwave::Vector();
        x_w = lightwave::Vector();
    }

    AffineVector(const lightwave::Vector &a, const lightwave::Vector &b,
                 const lightwave::Vector &c, const lightwave::Vector &d) {
        x_c = a;
        x_u = b;
        x_v = c;
        x_w = d;
    }

    AffineVector add(const AffineVector &other) {
        lightwave::Vector res_c = x_c + other.x_c;
        lightwave::Vector res_u = x_u + other.x_u;
        lightwave::Vector res_v = x_v + other.x_v;
        lightwave::Vector res_w = x_w + other.x_w;

        return AffineVector(res_c, res_u, res_v, res_w);
    }

    AffineVector dot(const AffineVector &other) {
        lightwave::Vector res_c = x_c * other.x_c;
        lightwave::Vector res_u = x_u * other.x_c + x_c * other.x_u;
        lightwave::Vector res_v = x_v * other.x_c + x_c * other.x_v;
        lightwave::Vector res_w =
            abs(x_w * other.x_c) + abs(x_c * other.x_w) +
            (abs(x_u) + abs(x_v) + x_w) *
                (abs(other.x_u) + abs(other.x_v) + other.x_w);

        return AffineVector(res_c, res_u, res_v, res_w);
    }

    AffineVector scalarMul(const lightwave::Vector4 &vec) {
        lightwave::Vector res_c = vec[0] * x_c;
        lightwave::Vector res_u = vec[1] * x_u;
        lightwave::Vector res_v = vec[2] * x_v;
        lightwave::Vector res_w = vec[3] * x_w;

        return AffineVector(res_c, res_u, res_v, res_w);
    }

    AffineVector matrixMul(const lightwave::Matrix3x3 &mat) {
        float a = mat(0, 0);
        float b = mat(0, 1);
        float c = mat(0, 2);

        float d = mat(1, 0);
        float e = mat(1, 1);
        float f = mat(1, 2);

        float g = mat(2, 0);
        float h = mat(2, 1);
        float i = mat(2, 2);

        lightwave::Vector res_c(a * x_c[0] + b * x_c[1] + c * x_c[2],
                                d * x_c[0] + e * x_c[1] + f * x_c[2],
                                g * x_c[0] + h * x_c[1] + i * x_c[2]);
        lightwave::Vector res_u(a * x_u[0] + b * x_u[1] + c * x_u[2],
                                d * x_u[0] + e * x_u[1] + f * x_u[2],
                                g * x_u[0] + h * x_u[1] + i * x_u[2]);
        lightwave::Vector res_v(a * x_v[0] + b * x_v[1] + c * x_v[2],
                                d * x_v[0] + e * x_v[1] + f * x_v[2],
                                g * x_v[0] + h * x_v[1] + i * x_v[2]);
        lightwave::Vector res_w(a * x_w[0] + b * x_w[1] + c * x_w[2],
                                d * x_w[0] + e * x_w[1] + f * x_w[2],
                                g * x_w[0] + h * x_w[1] + i * x_w[2]);

        return AffineVector(res_c, res_u, res_v, res_w);
    }
};

// MARK: D-BVH (Implicit Accelleration structure)

/// @brief A texel of the minmax-mipmap defined by 2 interger coordinates and a
/// mip level (LoD)
struct Texel {
    int LoD;
    Point2i ij;

    Texel() {
        LoD = 0;
        ij  = Point2i(0);
    }

    Texel(int startLoD, int i, int j) {
        LoD = startLoD;
        ij  = Point2i(i, j);
    }

    bool operator==(const Texel &other) const {
        return (LoD == other.LoD) && (ij == other.ij);
    }
    bool operator!=(const Texel &other) const { return !(*this == other); }

    /// @brief Goes one level down the mipmap chain
    void down() {
        LoD--;
        ij[0] = ij[0] * 2;
        ij[1] = ij[1] * 2;
    }

    /// @brief Goes one level up the mipmap chain
    void up() {
        LoD++;
        ij[0] = ij[0] / 2;
        ij[1] = ij[1] / 2;
    }

    /// @brief Moves to the next texel which needs to be tested for intersection
    /// or goes up if all texels at the current level are already tested
    void next() {
        while (1) {
            int i = 2 * (ij[0] % 2) + (ij[1] % 2);
            switch (i) {
            case 1:
                ij[0]++;
                ij[1]--;
                return;
            case 3:
                up();
                continue;
            default:
                ij[1]++;
                return;
            }
        }
    }
};

/// @brief  The minmax mipmap which acts as the implicit accelleration structure
struct MinMaxMipmap {
    std::vector<ref<Image>> m_images; // Images for the mipmaps. First color
                                      // channel is min and second color channel
                                      // is max
    int m_levels = 0;                 // Number of levels

    MinMaxMipmap() {};

    /// @brief Constructor generates the mipmap chain
    /// @param texture heightmap
    /// @param maxResolution the desired resolution of the heightmap
    /// @param cont
    MinMaxMipmap(ref<Texture> texture, const Point2i &maxResolution,
                 const Context &cont) {

        // Copy the minmax values for the first image
        m_images.push_back(std::make_shared<Image>(maxResolution));
        for (int y = 0; y < maxResolution.y(); y++) {
            for (int x = 0; x < maxResolution.x(); x++) {

                float a = texture->scalar(
                    Point2((x + 0.0f) / (float) maxResolution.x(),
                           (y + 0.0f) / (float) maxResolution.y()),
                    cont);
                float b = texture->scalar(
                    Point2((x + 1.0f) / (float) maxResolution.x(),
                           (y + 0.0f) / (float) maxResolution.y()),
                    cont);
                float c = texture->scalar(
                    Point2((x + 0.0f) / (float) maxResolution.x(),
                           (y + 1.0f) / (float) maxResolution.y()),
                    cont);
                float d = texture->scalar(
                    Point2((x + 1.0f) / (float) maxResolution.x(),
                           (y + 1.0f) / (float) maxResolution.y()),
                    cont);

                float minHeight = min(a, min(b, min(c, d)));
                float maxHeight = max(a, max(b, max(c, d)));

                m_images[0]->get(Point2i(x, y)) =
                    Color(minHeight, maxHeight, 0.0f);
            }
        }

        // Compute the desired number of levels and then compute them using
        // downsampling
        m_levels = static_cast<uint32_t>(std::floor(
                       std::log2(std::max(m_images[0]->resolution().x(),
                                          m_images[0]->resolution().y())))) +
                   1;
        int mipWidth  = m_images[0]->resolution().x();
        int mipHeight = m_images[0]->resolution().y();
        for (int i = 1; i < m_levels; i++) {
            m_images.push_back(std::make_shared<Image>(
                Point2i(mipWidth > 1 ? mipWidth / 2 : 1,
                        mipHeight > 1 ? mipHeight / 2 : 1)));

            for (int y = 0; y < (mipHeight > 1 ? mipHeight / 2 : 1); y++) {
                for (int x = 0; x < (mipWidth > 1 ? mipWidth / 2 : 1); x++) {
                    int x1 = (x * 2) + 0;
                    int y1 = (y * 2) + 0;
                    int x2 = x1 + 1 >= m_images[i - 1]->resolution().x()
                                 ? x1
                                 : x1 + 1;
                    int y2 = y1 + 1 >= m_images[i - 1]->resolution().y()
                                 ? y1
                                 : y1 + 1;

                    Color a = m_images[i - 1]->get(Point2i(x1, y1));
                    Color b = m_images[i - 1]->get(Point2i(x2, y1));
                    Color c = m_images[i - 1]->get(Point2i(x1, y2));
                    Color d = m_images[i - 1]->get(Point2i(x2, y2));

                    float minValue = min(a.r(), min(b.r(), min(c.r(), d.r())));
                    float maxValue = max(a.g(), max(b.g(), max(c.g(), d.g())));

                    float tl = texture->scalar(
                        Point2((x + 0.0f) / (float) (mipWidth / 2),
                               (y + 0.0f) / (float) (mipHeight / 2)),
                        cont,
                        i);
                    float tr = texture->scalar(
                        Point2((x + 1.0f) / (float) (mipWidth / 2),
                               (y + 0.0f) / (float) (mipHeight / 2)),
                        cont,
                        i);
                    float bl = texture->scalar(
                        Point2((x + 0.0f) / (float) (mipWidth / 2),
                               (y + 1.0f) / (float) (mipHeight / 2)),
                        cont,
                        i);
                    float br = texture->scalar(
                        Point2((x + 1.0f) / (float) (mipWidth / 2),
                               (y + 1.0f) / (float) (mipHeight / 2)),
                        cont,
                        i);

                    float minHeight = min(tl, min(tr, min(bl, br)));
                    float maxHeight = max(tl, max(tr, max(bl, br)));

                    minValue = min(minValue, minHeight);
                    maxValue = max(maxValue, maxHeight);

                    m_images[i]->get(Point2i(x, y)) =
                        Color(minValue, maxValue, 0.0f);
                }
            }

            mipWidth /= 2;
            mipHeight /= 2;
        }
    }
};

// MARK: - Displaced Mesh

/**
 * @brief A shape with a displacement map
 */
class DisplacedMesh : public AccelerationStructure {
    enum class IntersectionMode {
        Box,
        Triangulation,
    };

    /**
     * @brief The index buffer of the triangles.
     * The n-th element corresponds to the n-th triangle, and each component of
     * the element corresponds to one vertex index (into @c m_vertices ) of the
     * triangle. This list will always contain as many elements as there are
     * triangles.
     */
    std::vector<bob::TriangleIndices> m_triangles;
    /**
     * @brief The vertex buffer of the triangles, indexed by m_triangles.
     * Note that multiple triangles can share vertices, hence there can also be
     * fewer than @code 3 * numTriangles @endcode vertices.
     */
    std::vector<bob::Vertex> m_vertices;
    /// @brief The file this mesh was loaded from, for logging and debugging
    /// purposes.
    std::filesystem::path m_originalPath;
    /// @brief Whether to interpolate the normals from m_vertices, or report the
    /// geometric normal instead.
    bool m_smoothNormals;

    // Displacement Parameters
    float m_offset, m_scaling, m_bias;
    ref<Texture> m_displacementMap;

    MinMaxMipmap m_minmaxMipmap; // D-BVH
    int m_targetLoD;
    IntersectionMode m_intersectionMode;

    /// @brief We cache the inverse matrices and base mesh bounding boxes
    std::vector<Matrix3x3> m_invUVs;
    std::vector<Bounds> m_boudingBoxes;

protected:
    /// @brief Simple ray triangle intersection
    bool triangleIntersect(const Ray &ray, const Vertex &v0, const Vertex &v1,
                           const Vertex &v2, Intersection &its) const {
        // PROFILE("Triangle intersect");

        const auto edge1 = v1.position - v0.position;
        const auto edge2 = v2.position - v0.position;

        const auto pvec = ray.direction.cross(edge2);
        const auto det  = edge1.dot(pvec);
        if (det > -1e-8f && det < 1e-8f)
            return false;
        const auto invDet = 1 / det;

        const auto tvec = ray.origin - v0.position;
        const float u   = tvec.dot(pvec) * invDet;
        if (u < 0 || u > 1)
            return false;

        const auto qvec = tvec.cross(edge1);
        const float v   = ray.direction.dot(qvec) * invDet;
        if (v < 0 || u + v > 1)
            return false;

        const float t = edge2.dot(qvec) * invDet;
        if (t < Epsilon || t >= its.t)
            return false;

        const auto vInterpolated = Vertex::interpolate({ u, v }, v0, v1, v2);

        const auto deltaUV1 = v1.uv - v0.uv;
        const auto deltaUV2 = v2.uv - v0.uv;

        const float r =
            deltaUV1.x() * deltaUV2.y() - deltaUV1.y() * deltaUV2.x();
        const auto tangent =
            abs(r) < 1e-6 ? edge1
                          : (edge1 * deltaUV2.y() - edge2 * deltaUV1.y()) / r;

        its.t              = t;
        its.position       = vInterpolated.position;
        its.uv             = vInterpolated.uv;
        its.geometryNormal = edge1.cross(edge2).normalized();
        its.shadingNormal  = m_smoothNormals ? vInterpolated.normal.normalized()
                                             : its.geometryNormal;
        its.tangent        = tangent;
        its.pdf            = 0;

        assert_finite(its.t, {
            logger(EError,
                   "triangle: %s %s %s",
                   v0.position,
                   v1.position,
                   v2.position);
            logger(EError, "count: %d", m_vertices.size());
            logger(EError, "offending shape: %s", this);
        });
        return true;
    }

    int numberOfPrimitives() const override { return int(m_triangles.size()); }

    /// @brief Test for ray - bounding box intersection using the slab method
    bool intersectAABB(const Bounds &aabb, const Ray &ray,
                       Intersection &its) const {
        // PROFILE("AABB intersect");

        // Slab method
        float t1, t2, t_near = -FLT_MAX, t_far = FLT_MAX;
        int hitAxis = -1;

        for (int i = 0; i < 3; ++i) {
            if (ray.direction[i] != 0.0f) {
                t1 = (aabb.min()[i] - ray.origin[i]) / ray.direction[i];
                t2 = (aabb.max()[i] - ray.origin[i]) / ray.direction[i];

                if (t1 > t2) {
                    float temp = t1;
                    t1         = t2;
                    t2         = temp;
                }

                if (t1 > t_near) {
                    t_near  = t1;
                    hitAxis = i;
                }
                if (t2 < t_far)
                    t_far = t2;

                if (t_near > t_far || t_far < 0.0f) {
                    return false; // No intersection
                }
            } else {
                if (ray.origin[i] < aabb.min()[i] ||
                    ray.origin[i] > aabb.max()[i]) {
                    return false; // No intersectio
                }
            }
        }

        Vector intNormal(0.0f);
        Vector tangent(0.0f);
        if (hitAxis != -1) {
            intNormal[0] = intNormal[1] = intNormal[2] = 0.0f;
            if (fabs(t_near - (aabb.min()[hitAxis] - ray.origin[hitAxis]) /
                                  ray.direction[hitAxis]) < 1e-6) {
                intNormal[hitAxis] = -1.0f;
            } else {
                intNormal[hitAxis] = 1.0f;
            }

            Vector up = Vector(0.0f, 1.0f, 0.0f);
            if (intNormal[1] == 1.0f || intNormal[1] == -1.0f) {
                up[0] =
                    1.0f; // Use a different up vector if normal is along Y axis
                up[1] = 0.0f;
            }

            tangent = intNormal.cross(up);
            tangent = tangent.normalized();
        }

        float t = t_near;

        if (t < Epsilon || t >= its.t)
            return false;

        its.t              = t;
        its.position       = ray(t);
        its.uv             = Point2();
        its.geometryNormal = intNormal.normalized();
        its.shadingNormal  = its.geometryNormal;
        its.tangent        = tangent;
        its.pdf            = 0;

        return true;
    }

    /// @brief Uses the algorithm in the paper and affine arithmetic to compute
    /// an Bounding Box on the fly for a given texel of the minmax-Mipmap
    Bounds buildAABBFromTexel(int primitiveIndex, const Texel &texel) const {
        // PROFILE("AABB construction");
        const auto &triangle = m_triangles[primitiveIndex];
        const auto &vertex0  = m_vertices[triangle[0]];
        const auto &vertex1  = m_vertices[triangle[1]];
        const auto &vertex2  = m_vertices[triangle[2]];

        lightwave::Vector p1(vertex0.position.data()),
            p2(vertex1.position.data()), p3(vertex2.position.data());
        lightwave::Vector n1(vertex0.normal.decompress().data()),
            n2(vertex1.normal.decompress().data()),
            n3(vertex2.normal.decompress().data());

        lightwave::Matrix3x3 t_uv = m_invUVs[primitiveIndex];
        lightwave::Matrix3x3 pos(
            { p1[0], p2[0], p3[0], p1[1], p2[1], p3[1], p1[2], p2[2], p3[2] });
        lightwave::Matrix3x3 normal(
            { n1[0], n2[0], n3[0], n1[1], n2[1], n3[1], n1[2], n2[2], n3[2] });

        float W  = m_minmaxMipmap.m_images[0]->resolution().x(),
              H  = m_minmaxMipmap.m_images[0]->resolution().y();
        float ii = texel.ij[0], j = texel.ij[1], k = texel.LoD;
        lightwave::Vector x_c =
            powf(2.0f, k) *
            lightwave::Vector((ii + 0.5f) / W, (j + 0.5f) / H, 1.0f);
        x_c[2] = 1.0f;
        lightwave::Vector x_u =
            powf(2.0f, k) * lightwave::Vector(1.0f / (2.0f * W), 0.0f, 0.0f);
        lightwave::Vector x_v =
            powf(2.0f, k) * lightwave::Vector(0.0f, 1.0f / (2.0f * H), 0.0f);
        lightwave::Vector x_w =
            powf(2.0f, k) * lightwave::Vector(0.0f, 0.0f, 0.0f);

        lightwave::Matrix3x3 posTransform    = pos * t_uv;
        lightwave::Matrix3x3 normalTransform = normal * t_uv;

        Color minmax = m_minmaxMipmap.m_images[k]->get(Point2i(ii, j));

        float min = m_offset + m_scaling * (minmax.r() - m_bias);
        float max = m_offset + m_scaling * (minmax.g() - m_bias);
        lightwave::Vector h_c =
            0.5f * lightwave::Vector(max + min, max + min, max + min);
        lightwave::Vector h_u = lightwave::Vector(0, 0, 0);
        lightwave::Vector h_v = lightwave::Vector(0, 0, 0);
        lightwave::Vector h_w =
            0.5f * lightwave::Vector(max - min, max - min, max - min);
        AffineVector h(h_c, h_u, h_v, h_w);

        AffineVector uv(x_c, x_u, x_v, x_w);

        AffineVector affinePos    = uv.matrixMul(posTransform);
        AffineVector affineNormal = uv.matrixMul(normalTransform);

        AffineVector displacedSurface = affinePos.add(h.dot(affineNormal));

        Bounds aabb(
            Point(displacedSurface.x_c - abs(displacedSurface.x_u) -
                  abs(displacedSurface.x_v) - abs(displacedSurface.x_w)),
            Point(displacedSurface.x_c + abs(displacedSurface.x_u) +
                  abs(displacedSurface.x_v) + abs(displacedSurface.x_w)));

        return aabb;
    }

    /// @brief Custom intersect function which uses the algorithm from the paper
    /// to account for displacements
    bool intersect(int primitiveIndex, const Ray &ray, Intersection &its,
                   Sampler &rng, Context &cont) const override {
        PROFILE("Base mesh intersect");

        const auto &triangle = m_triangles[primitiveIndex];
        const auto &vertex0  = m_vertices[triangle[0]];
        const auto &vertex1  = m_vertices[triangle[1]];
        const auto &vertex2  = m_vertices[triangle[2]];

        Matrix2x2 r = Matrix2x2({ 0.0f, -1.0f, 1.0f, 0.0f });
        Vector2 p1(vertex0.uv.data());
        Vector2 p2(vertex1.uv.data());
        Vector2 p3(vertex2.uv.data());
        Vector2 n1 = (r * (p1 - p2)).normalized();
        Vector2 n2 = (r * (p2 - p3)).normalized();
        Vector2 n3 = (r * (p3 - p1)).normalized();

        bool intersectionFound = false;

        // Start at root
        Texel texel    = Texel(m_minmaxMipmap.m_levels - 1, 0, 0);
        Texel endTexel = Texel(m_minmaxMipmap.m_levels - 1, 0, 0);
        endTexel.next();

        Intersection hitAABBIts = Intersection();
        Texel hitTexel;
        Vector2 hitTexelCenter;
        float hitTexelHalfWidth;

        // Walk through the D-BVH (minmax mipmap), create bounding boxes and
        // intersect with them
        while (texel != endTexel) {
            float texelHalfWidth =
                (1.0f / m_minmaxMipmap.m_images[texel.LoD]->resolution().x()) *
                0.5f;
            Vector2 texelCenter =
                Vector2((float) texel.ij.x() + 0.5f,
                        (float) texel.ij.y() + 0.5f) /
                Vector2(m_minmaxMipmap.m_images[texel.LoD]->resolution().x());
            CollisionType collisionResult =
                Collision(p1, p2, p3, n1, n2, n3, texelCenter, texelHalfWidth);

            if (collisionResult == SquareOutsideTriangle) {
                texel.next();
            } else {
                Bounds aabb = buildAABBFromTexel(primitiveIndex, texel);
                Intersection tempIntersection;
                if (intersectAABB(aabb, ray, tempIntersection) == false) {
                    texel.next();
                } else if (texel.LoD <= m_targetLoD) {

                    // We are at an leaf node so we have to choose an
                    // intersection method

                    if (tempIntersection.t < hitAABBIts.t) {
                        // Box intersection
                        hitAABBIts = tempIntersection;
                        intersectionFound |=
                            (m_intersectionMode == IntersectionMode::Box);
                    }

                    // With intersection mode triangulation we create two temp
                    // triangles
                    if (m_intersectionMode == IntersectionMode::Triangulation) {
                        hitTexel          = texel;
                        hitTexelCenter    = texelCenter;
                        hitTexelHalfWidth = texelHalfWidth;

                        auto heightImage =
                            m_minmaxMipmap.m_images[hitTexel.LoD];
                        Point2 hitMipmapSize =
                            heightImage->resolution().cast<float>();
                        Point2 texelCoordsTL =
                            hitTexel.ij.cast<float>() + Vector2(0.0f, 0.0f);
                        Point2 texelCoordsTR =
                            hitTexel.ij.cast<float>() + Vector2(1.0f, 0.0f);
                        Point2 texelCoordsBL =
                            hitTexel.ij.cast<float>() + Vector2(0.0f, 1.0f);
                        Point2 texelCoordsBR =
                            hitTexel.ij.cast<float>() + Vector2(1.0f, 1.0f);
                        float heightTL = m_displacementMap->scalar(
                            Point2(texelCoordsTL.x() / hitMipmapSize.x(),
                                   texelCoordsTL.y() / hitMipmapSize.y()),
                            cont,
                            hitTexel.LoD);
                        float heightTR = m_displacementMap->scalar(
                            Point2(texelCoordsTR.x() / hitMipmapSize.x(),
                                   texelCoordsTR.y() / hitMipmapSize.y()),
                            cont,
                            hitTexel.LoD);
                        float heightBL = m_displacementMap->scalar(
                            Point2(texelCoordsBL.x() / hitMipmapSize.x(),
                                   texelCoordsBL.y() / hitMipmapSize.y()),
                            cont,
                            hitTexel.LoD);
                        float heightBR = m_displacementMap->scalar(
                            Point2(texelCoordsBR.x() / hitMipmapSize.x(),
                                   texelCoordsBR.y() / hitMipmapSize.y()),
                            cont,
                            hitTexel.LoD);

                        Vector texelTL =
                            Vector(hitTexelCenter.x() - hitTexelHalfWidth,
                                   hitTexelCenter.y() - hitTexelHalfWidth,
                                   1.0f);
                        Vector texelTR =
                            Vector(hitTexelCenter.x() + hitTexelHalfWidth,
                                   hitTexelCenter.y() - hitTexelHalfWidth,
                                   1.0f);
                        Vector texelBL =
                            Vector(hitTexelCenter.x() - hitTexelHalfWidth,
                                   hitTexelCenter.y() + hitTexelHalfWidth,
                                   1.0f);
                        Vector texelBR =
                            Vector(hitTexelCenter.x() + hitTexelHalfWidth,
                                   hitTexelCenter.y() + hitTexelHalfWidth,
                                   1.0f);

                        float u1 = vertex0.uv[0], u2 = vertex1.uv[0],
                              u3 = vertex2.uv[0];
                        float v1 = vertex0.uv[1], v2 = vertex1.uv[1],
                              v3 = vertex2.uv[1];

                        lightwave::Vector pp1(vertex0.position.data()),
                            pp2(vertex1.position.data()),
                            pp3(vertex2.position.data());
                        lightwave::Vector nn1(
                            vertex0.normal.decompress().data()),
                            nn2(vertex1.normal.decompress().data()),
                            nn3(vertex2.normal.decompress().data());

                        lightwave::Matrix3x3 t_uv_inv(
                            { u1, u2, u3, v1, v2, v3, 1, 1, 1 });
                        lightwave::Matrix3x3 t_uv =
                            invert(t_uv_inv).value_or(Matrix3x3::identity());
                        lightwave::Matrix3x3 pos({ pp1[0],
                                                   pp2[0],
                                                   pp3[0],
                                                   pp1[1],
                                                   pp2[1],
                                                   pp3[1],
                                                   pp1[2],
                                                   pp2[2],
                                                   pp3[2] });
                        lightwave::Matrix3x3 normal({ nn1[0],
                                                      nn2[0],
                                                      nn3[0],
                                                      nn1[1],
                                                      nn2[1],
                                                      nn3[1],
                                                      nn1[2],
                                                      nn2[2],
                                                      nn3[2] });

                        Vector triangleNormal1 = normal * t_uv * texelTL;
                        Vector triangleNormal2 = normal * t_uv * texelTR;
                        Vector triangleNormal3 = normal * t_uv * texelBL;
                        Vector triangleNormal4 = normal * t_uv * texelBR;

                        Vector trianglePos1 =
                            (pos * t_uv * texelTL) +
                            (m_offset + m_scaling * heightTL) * triangleNormal1;
                        Vector trianglePos2 =
                            (pos * t_uv * texelTR) +
                            (m_offset + m_scaling * heightTR) * triangleNormal2;
                        Vector trianglePos3 =
                            (pos * t_uv * texelBL) +
                            (m_offset + m_scaling * heightBL) * triangleNormal3;
                        Vector trianglePos4 =
                            (pos * t_uv * texelBR) +
                            (m_offset + m_scaling * heightBR) * triangleNormal4;

                        Vertex vert0, vert1, vert2, vert3;
                        vert0.position = Point(trianglePos1.x(),
                                               trianglePos1.y(),
                                               trianglePos1.z());
                        vert1.position = Point(trianglePos2.x(),
                                               trianglePos2.y(),
                                               trianglePos2.z());
                        vert2.position = Point(trianglePos3.x(),
                                               trianglePos3.y(),
                                               trianglePos3.z());
                        vert3.position = Point(trianglePos4.x(),
                                               trianglePos4.y(),
                                               trianglePos4.z());
                        vert0.normal   = triangleNormal1.normalized();
                        vert1.normal   = triangleNormal2.normalized();
                        vert2.normal   = triangleNormal3.normalized();
                        vert3.normal   = triangleNormal4.normalized();
                        vert0.uv       = Vector2(texelTL.x(), texelTL.y());
                        vert1.uv       = Vector2(texelTR.x(), texelTR.y());
                        vert2.uv       = Vector2(texelBL.x(), texelBL.y());
                        vert3.uv       = Vector2(texelBR.x(), texelBR.y());

                        intersectionFound |=
                            triangleIntersect(ray, vert0, vert1, vert2, its);
                        intersectionFound |=
                            triangleIntersect(ray, vert1, vert3, vert2, its);
                        its.shadingNormal = its.geometryNormal;
                    }

                    texel.next();
                } else {
                    texel.down();
                }
            }
        }

        if (!intersectionFound) {
            return false;
        }

        if (m_intersectionMode == IntersectionMode::Box) {
            its.t              = hitAABBIts.t;
            its.position       = ray(hitAABBIts.t);
            its.uv             = Point2();
            its.geometryNormal = hitAABBIts.geometryNormal;
            its.shadingNormal  = hitAABBIts.shadingNormal;
            its.tangent        = hitAABBIts.tangent;
            its.pdf            = 0;
        }

        return true;
    }

    Bounds getBoundingBox(int primitiveIndex) const override {
        return m_boudingBoxes[primitiveIndex];
    }

    Point getCentroid(int primitiveIndex) const override {
        return getBoundingBox(primitiveIndex).center();
    }

    void loadMesh(const std::filesystem::path path) {
        bob::Mesh mesh(path);
        m_vertices  = mesh.getVertices();
        m_triangles = mesh.getTriangles();
    }

    void reorderPrimitives(const std::vector<int> &newOrder) override {
        std::vector<bob::TriangleIndices> updatedTriangles;
        for (int index : newOrder) {
            updatedTriangles.push_back(m_triangles[index]);
        }
        m_triangles = updatedTriangles;
    }

public:
    DisplacedMesh(const Properties &properties) {
        m_originalPath  = properties.get<std::filesystem::path>("filename");
        m_smoothNormals = properties.get<bool>("smooth", true);
        if (m_originalPath.extension() == ".bob")
            loadMesh(m_originalPath);
        else
            readPLY(m_originalPath, m_triangles, m_vertices);
        logger(EInfo,
               "loaded ply with %d triangles, %d vertices",
               m_triangles.size(),
               m_vertices.size());

        // Read the displacement paramters
        Context cont;
        m_displacementMap = properties.get<Texture>("height");
        m_minmaxMipmap =
            MinMaxMipmap(m_displacementMap, Point2i(2048, 2048), cont);
        m_offset    = properties.get<Texture>("offset")->scalar(Point2(), cont);
        m_scaling   = properties.get<Texture>("scale")->scalar(Point2(), cont);
        m_bias      = properties.get<float>("bias", 0.0f);
        m_targetLoD = properties.get<int>("lod", 0);
        m_intersectionMode = properties.getEnum<IntersectionMode>(
            "intersection",
            IntersectionMode::Triangulation,
            {
                { "box", IntersectionMode::Box },
                { "triangulation", IntersectionMode::Triangulation },
            });

        // Precompute inverse matrices and the base mesh bounding boxes
        for (int primitiveIndex = 0; primitiveIndex < (int) m_triangles.size();
             primitiveIndex++) {
            const auto &triangle = m_triangles[primitiveIndex];
            const auto &vertex0  = m_vertices[triangle[0]];
            const auto &vertex1  = m_vertices[triangle[1]];
            const auto &vertex2  = m_vertices[triangle[2]];

            float u1 = vertex0.uv[0], u2 = vertex1.uv[0], u3 = vertex2.uv[0];
            float v1 = vertex0.uv[1], v2 = vertex1.uv[1], v3 = vertex2.uv[1];

            lightwave::Matrix3x3 t_uv({ u1, u2, u3, v1, v2, v3, 1, 1, 1 });
            m_invUVs.push_back(invert(t_uv).value_or(Matrix3x3::identity()));

            Matrix2x2 r = Matrix2x2({ 0.0f, -1.0f, 1.0f, 0.0f });
            Vector2 p1(vertex0.uv.data());
            Vector2 p2(vertex1.uv.data());
            Vector2 p3(vertex2.uv.data());
            Vector2 n1 = (r * (p1 - p2)).normalized();
            Vector2 n2 = (r * (p2 - p3)).normalized();
            Vector2 n3 = (r * (p3 - p1)).normalized();

            // Start at root
            Texel texel    = Texel(m_minmaxMipmap.m_levels - 1, 0, 0);
            Texel endTexel = Texel(m_minmaxMipmap.m_levels - 1, 0, 0);
            endTexel.next();

            Bounds bounds = Bounds();

            while (texel != endTexel) {
                float texelHalfWidth =
                    (1.0f /
                     m_minmaxMipmap.m_images[texel.LoD]->resolution().x()) *
                    0.5f;
                Vector2 texelCenter =
                    Vector2((float) texel.ij.x() + 0.5f,
                            (float) texel.ij.y() + 0.5f) /
                    Vector2(
                        m_minmaxMipmap.m_images[texel.LoD]->resolution().x());
                CollisionType collisionResult = Collision(
                    p1, p2, p3, n1, n2, n3, texelCenter, texelHalfWidth);

                if (collisionResult == SquareOutsideTriangle) {
                    texel.next();
                } else if (collisionResult == SquareInsideTriangle ||
                           texel.LoD == 0) {
                    Bounds aabb = buildAABBFromTexel(primitiveIndex, texel);

                    bounds.extend(aabb);
                    texel.next();
                } else {
                    texel.down();
                }
            }

            m_boudingBoxes.push_back(bounds);
        }

        buildAccelerationStructure();
    }

    bool intersect(const Ray &ray, Intersection &its, Sampler &rng,
                   Context &cont) const override {
        PROFILE("Displaced mesh");
        return AccelerationStructure::intersect(ray, its, rng, cont);
    }

    AreaSample sampleArea(Sampler &rng, Context &cont) const override{
        // only implement this if you need triangle mesh area light sampling for
        // your rendering competition
        NOT_IMPLEMENTED
    }

    std::string toString() const override {
        return tfm::format(
            "DisplacedMesh[\n"
            "  vertices = %d,\n"
            "  triangles = %d,\n"
            "  filename = \"%s\"\n"
            "]",
            m_vertices.size(),
            m_triangles.size(),
            m_originalPath.generic_string());
    }
};

} // namespace lightwave

REGISTER_SHAPE(DisplacedMesh, "displaced")
