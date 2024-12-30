#pragma once

#include "math.hpp"
#include "quantization.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace internal {

class TriangleMesh {

    std::vector<Vector3i> m_triangles;
    std::vector<Vertex> m_vertices;

    bool has_uvs;

private:
    void generate_uvs() {
        Bounds bbox;
        for (const Vertex &v : m_vertices)
            bbox.extend(v.position);

        for (size_t i = 0; i < m_vertices.size(); ++i) {
            auto &v        = m_vertices.at(i);
            const Vector d = bbox.diagonal();
            const Vector t = v.position - bbox.min();

            Vector2 p = Vector2(0);
            if (d.x() > Epsilon)
                p.x() = t.x() / d.x();
            if (d.y() > Epsilon)
                p.y() = t.y() / d.y();
            v.uv = p; // Drop the z coordinate
        }
    }

public:
    TriangleMesh() = default;
    TriangleMesh(const bool uvs) : has_uvs(uvs) {}

    void addVertex(const Vertex vs) { m_vertices.push_back(vs); }
    void addFace(const Vector3i indices) { m_triangles.push_back(indices); }

    void saveMesh(const std::string path) {
        if (!has_uvs)
            generate_uvs();

        std::ofstream file;

        // clear contents of file if it already exists
        file.open(path, std::ios::out | std::ios::binary);
        file.close();

        // append data to empty file
        file.open(path, std::ios::app | std::ios::binary);

        size_t numVertices  = m_vertices.size();
        size_t numTriangles = m_triangles.size();

        file.write((char *) &numVertices, sizeof(numVertices));
        file.write((char *) &numTriangles, sizeof(numTriangles));

        for (Vertex v : m_vertices) {
            v.write(file);
        }

        file.write((char *) &m_triangles[0],
                   m_triangles.size() * sizeof(Vector3i));

        file.close();
    }
};

} // namespace internal