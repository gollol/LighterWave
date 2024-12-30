#include "math.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

enum ObjectType { UNKNOWN = -1, FUR = 0, MESH = 1 };

namespace internal {

class Hair {
protected:
    std::vector<Vector> m_keyPoints;

public:
    void addPoint(const Vector &keyPoint) { m_keyPoints.push_back(keyPoint); }
    std::vector<Vector> toBezier() const {
        int numCp = m_keyPoints.size();
        std::vector<Vector> temp(numCp + 2);

        temp[0]         = m_keyPoints[0];
        temp[numCp + 1] = m_keyPoints[numCp - 1];

        for (int i = 0; i < numCp; i++)
            temp[i + 1] = m_keyPoints[i];

        std::vector<Vector> bezier(3 * numCp - 2);
        bezier[0] = temp[1];
        for (int i = 0; i < numCp - 1; i++) {
            bezier[3 * i + 1] = temp[i + 1] + (temp[i + 2] - temp[i]) / 6;
            bezier[3 * i + 2] = temp[i + 2] - (temp[i + 3] - temp[i + 1]) / 6;
            bezier[3 * i + 3] = temp[i + 2];
        }
        return bezier;
    }
};

class Fur {
protected:
    std::vector<Vector> m_keyPoints;
    int m_curvePoints;

public:
    Fur() = default;
    Fur(int keyPoints) : m_curvePoints(3 * keyPoints - 2) {}

    void addHair(const Hair &h) {
        std::vector<Vector> bezierPoints = h.toBezier();
        m_keyPoints.insert(
            m_keyPoints.end(), bezierPoints.begin(), bezierPoints.end());
    }

    void saveFur(const std::string path) {
        std::ofstream file;

        // clear contents of file if it already exists
        file.open(path, std::ios::out | std::ios::binary);
        file.close();

        // append data to empty file
        file.open(path, std::ios::app | std::ios::binary);

        ObjectType type = ObjectType::FUR;
        int size        = m_keyPoints.size();
        file.write((char *) &type, sizeof(ObjectType));
        file.write((char *) &size, sizeof(size));

        file.write((char *) &m_curvePoints, sizeof(m_curvePoints));

        file.write((char *) &m_keyPoints[0],
                   m_keyPoints.size() * sizeof(Vector));

        file.close();
    }
};
} // namespace internal