#include <lightwave/math.hpp>

#include <bob.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace lightwave {

void readPLY(const std::filesystem::path &path,
             std::vector<bob::TriangleIndices> &indices,
             std::vector<bob::Vertex> &vertices);

} // namespace lightwave
