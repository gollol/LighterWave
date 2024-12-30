#pragma once
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "src/math.hpp"
#include "src/quantization.hpp"
#include "src/utils.hpp"

namespace bob {

class Vector3 : protected internal::Vector {
public:
    Vector3() = default;
    Vector3(const float x, const float y, const float z)
        : internal::Vector(x, y, z) {}
    Vector3(internal::Vector v) : internal::Vector(v) {}
    /// @brief returns the components of the vector,
    /// sequentially storred in an array
    std::array<float, 3> data() const { return m_data; }
    const float &operator[](int i) const { return m_data[i]; }
    float &operator[](int i) { return m_data[i]; }
};

class Vector3i : protected internal::Vector3i {
public:
    Vector3i() = default;
    Vector3i(const int x, const int y, const int z)
        : internal::Vector3i(x, y, z) {}
    /// @brief returns the components of the vector,
    /// sequentially storred in an array
    std::array<int, 3> data() const { return m_data; }
    const int &operator[](int i) const { return m_data[i]; }
    int &operator[](int i) { return m_data[i]; }
};

class Vector2 : protected internal::Vector2 {
public:
    Vector2() = default;
    Vector2(const float x, const float y) : internal::Vector2(x, y) {}
    Vector2(std::array<float, 2> v) : internal::Vector2(v) {}
    /// @brief returns the components of the vector,
    /// sequentially storred in an array
    std::array<float, 2> data() const { return m_data; }
    const float &operator[](int i) const { return m_data[i]; }
    float &operator[](int i) { return m_data[i]; }
};

class CompressedNormal : protected internal::Oct32 {
public:
    CompressedNormal() = default;
    CompressedNormal(const internal::Oct32 &v) : internal::Oct32(v) {}
    CompressedNormal(const std::array<float, 3> &v)
        : internal::Oct32(internal::Vector(v)) {}
    /// @brief returns the components of the vector,
    /// sequentially storred in an array
    Vector3 decompress() const { return Vector3(internal::Oct32(m_encoded)); }
    std::array<int16_t, 2> data() const {
        return { m_encoded.data()[0].data(), m_encoded.data()[1].data() };
    }
};

typedef Vector3 Point;
typedef Vector3i TriangleIndices;

struct Vertex {
    Point position;
    Vector2 uv;
    CompressedNormal normal;

    Vertex() = default;
    Vertex(internal::Vertex v)
        : position(v.position.data()), uv(v.uv.data()), normal(v.normal) {}

    Vertex(std::ifstream &file) {
        file.read((char *) &position, sizeof(position));
        file.read((char *) &uv, sizeof(uv));
        file.read((char *) &normal, sizeof(normal));
    }
};

class Fur {
    std::vector<Vector3> m_keyPoints;
    int m_curvePoints;

public:
    /// @brief loads a .bob file that contains a Fur object from disk
    Fur(const std::filesystem::path &path) {
        std::ifstream file;
        file.open(path, std::ios::in | std::ios::binary);

        ObjectType type;
        int size;
        file.read((char *) &type, sizeof(ObjectType));
        if (type != ObjectType::FUR) {
            std::cerr << "File is not a Fur file." << std::endl;
            return;
        }
        file.read((char *) &size, sizeof(size));

        file.read((char *) &m_curvePoints, sizeof(m_curvePoints));

        m_keyPoints.resize(size);
        file.read((char *) &m_keyPoints[0], sizeof(Vector3) * size);

        file.close();
    }

    /// @brief returns the key points of all hair sequentially stored in a
    /// vector
    std::vector<Vector3> &getKeyPoints() { return m_keyPoints; }
    /// @brief returns the number of key points per curve
    int getKeyPointCount() { return m_curvePoints; }
};

class Mesh {
    std::vector<Vertex> m_vertices;
    std::vector<TriangleIndices> m_triangles;

public:
    Mesh(const std::filesystem::path &path) {
        std::ifstream file;
        file.open(path, std::ios::in | std::ios::binary);

        size_t numVertices, numTriangles;
        file.read((char *) &numVertices, sizeof(numVertices));
        file.read((char *) &numTriangles, sizeof(numTriangles));

        m_vertices.reserve(numVertices);
        m_triangles.resize(numTriangles);

        for (size_t i = 0; i < numVertices; ++i) {
            m_vertices.emplace_back(file);
        }

        file.read((char *) &m_triangles[0], sizeof(Vector3i) * numTriangles);
        file.close();
    }

    std::vector<Vertex> getVertices() { return m_vertices; }
    std::vector<Vector3i> getTriangles() { return m_triangles; }
};
} // namespace bob