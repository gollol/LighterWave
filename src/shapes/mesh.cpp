#include <lightwave.hpp>

#include "../core/plyparser.hpp"
#include "accel.hpp"
#include "bob.hpp"
#include <fstream>

namespace lightwave {

/**
 * @brief A shape consisting of many (potentially millions) of triangles, which
 * share an index and vertex buffer. Since individual triangles are rarely
 * needed (and would pose an excessive amount of overhead), collections of
 * triangles are combined in a single shape.
 */
class TriangleMesh : public AccelerationStructure {
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

protected:
    int numberOfPrimitives() const override { return int(m_triangles.size()); }

    bool intersect(int primitiveIndex, const Ray &ray, Intersection &its,
                   Sampler &rng, Context &cont) const override {
        const auto &triangle = m_triangles[primitiveIndex];
        const Vertex v0      = m_vertices[triangle[0]];
        const Vertex v1      = m_vertices[triangle[1]];
        const Vertex v2      = m_vertices[triangle[2]];

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
        if (!its.testAlpha(vInterpolated.uv, rng, cont))
            return false;

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
        // its.shadingNormal = m_smoothNormals ?
        // decode_spherical(vInterpolated.normal).normalized() :
        // its.geometryNormal;
        its.shadingNormal = m_smoothNormals ? vInterpolated.normal.normalized()
                                            : its.geometryNormal;
        its.tangent       = tangent;
        its.pdf           = 0;

        cont.curvesIndex = -1;

        assert_finite(its.t, {
            logger(EError,
                   "triangle: %s %s %s",
                   v0.position,
                   v1.position,
                   v2.position);
            logger(EError,
                   "indices: %d %d %d",
                   triangle[0],
                   triangle[1],
                   triangle[2]);
            logger(EError, "count: %d", m_vertices.size());
            logger(EError, "offending shape: %s", this);
        });
        return true;
        // hints:
        // * use m_triangles[primitiveIndex] to get the vertex indices of the
        // triangle that should be intersected
        // * if m_smoothNormals is true, interpolate the vertex normals from
        // m_vertices
        //   * make sure that your shading frame stays orthonormal!
        // * if m_smoothNormals is false, use the geometrical normal (can be
        // computed from the vertex positions)
    }

    Bounds getBoundingBox(int primitiveIndex) const override {
        const auto &triangle = m_triangles[primitiveIndex];
        const auto &v0       = m_vertices[triangle[0]];
        const auto &v1       = m_vertices[triangle[1]];
        const auto &v2       = m_vertices[triangle[2]];

        Bounds result;
        result.extend({ v0.position.data() });
        result.extend({ v1.position.data() });
        result.extend({ v2.position.data() });
        return result;
    }

    Point getCentroid(int primitiveIndex) const override {
        const auto &triangle = m_triangles[primitiveIndex];
        const bob::Vertex v0 = m_vertices[triangle[0]];
        const bob::Vertex v1 = m_vertices[triangle[1]];
        const bob::Vertex v2 = m_vertices[triangle[2]];
        return (Vector(v0.position.data()) + Vector(v1.position.data()) +
                Vector(v2.position.data())) /
               3;
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
    TriangleMesh(const Properties &properties) {
        m_originalPath  = properties.get<std::filesystem::path>("filename");
        m_smoothNormals = properties.get<bool>("smooth", true);
        logger(EInfo, "start loading %s", m_originalPath.filename().string());

        if (m_originalPath.extension() == ".bob")
            loadMesh(m_originalPath);
        else
            readPLY(m_originalPath, m_triangles, m_vertices);

        logger(EInfo,
               "loaded %s file with %d triangles, %d vertices",
               m_originalPath.extension().string(),
               m_triangles.size(),
               m_vertices.size());
        buildAccelerationStructure();
    }

    bool intersect(const Ray &ray, Intersection &its, Sampler &rng,
                   Context &cont) const override {
        PROFILE("Triangle mesh")
        return AccelerationStructure::intersect(ray, its, rng, cont);
    }

    AreaSample sampleArea(Sampler &rng, Context &cont) const override{
        // only implement this if you need triangle mesh area light sampling for
        // your rendering competition
        NOT_IMPLEMENTED
    }

    std::string toString() const override {
        return tfm::format(
            "Mesh[\n"
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

REGISTER_SHAPE(TriangleMesh, "mesh")
