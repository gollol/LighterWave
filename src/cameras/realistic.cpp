#include <fstream>
#include <lightwave.hpp>

namespace fs = std::filesystem;
using namespace std::string_literals;
/**
 * @brief Implements a realistic camera model as described in the PBRT book:
 * https://www.pbr-book.org/3ed-2018/Camera_Models/Realistic_Cameras
 */
namespace lightwave {

class Realistic : public Camera {

public:
    /// @brief Height and width of the sensor (m).
    Vector2 m_sensorSize;

    /// @brief Desired depth of the plane of focus (m).
    float m_focusDistance;

    /// @brief The lenses that make up the camera.
    std::vector<ref<Lens>> m_elementInterfaces;

    Realistic(const Properties &properties) : Camera(properties) {
        m_sensorSize.x() = 0.001 * properties.get<float>("sensorWidth", 36);
        m_sensorSize.y() = 0.001 * properties.get<float>("sensorHeight", 24);
        m_focusDistance  = properties.get<float>("focusDistance");

        m_elementInterfaces = properties.getChildren<Lens>();

        // Draw lens system
        if (properties.has("pathToSchematic")) {
            auto filePath = properties.get<fs::path>("pathToSchematic");
            std::ofstream svgFile(filePath);
            svgFile << createHeader() << std::endl;
            svgFile << drawOpticalAxis() << std::endl;

            float zPos = -1000 * lensFrontZ() + 10;
            for (size_t i = 0; i < m_elementInterfaces.size(); i++) {
                svgFile << drawInterface(&zPos, i) << std::endl;
            }
            svgFile << "</svg>" << std::endl;
            logger(EInfo, "Save lens schematic to: %s", filePath);
            svgFile.close();
        }

        // Adjust distance between lens stack and film for given focus distance
        m_elementInterfaces.front()->setThickness(focusThickLens());
    }

    /// @brief Sample point on film plane and trace from it through camera.
    CameraSample sample(const Point2 &normalized, Sampler &rng) const override {

        Point filmPoint = Point(-m_sensorSize.x() * normalized.x(),
                                -m_sensorSize.y() * normalized.y(),
                                0);

        const float frontRadius = frontLensRadius();
        Point2 lensCoords       = squareToUniformDiskConcentric(rng.next2D());
        Point lensPoint         = Point(lensCoords.x() * frontRadius,
                                lensCoords.y() * frontRadius,
                                lensFrontZ());

        Vector dirFilmRay = (lensPoint - filmPoint).normalized();

        Ray rFilm = Ray(filmPoint, dirFilmRay).normalized();
        Ray outRay;

        if (!traceLensesFromFilm(rFilm, outRay))
            return CameraSample::invalid();

        const float time = rng.next();
        outRay.time = time;
        SimpleTransform T = m_transform->interpolate(time);
        
        return CameraSample{
            .ray    = T.apply(outRay).normalized(),
            .weight = Color(1.f),
        };
    }

    std::string toString() const override {
        return tfm::format(
            "Realistic[\n"
            "  width = %d,\n"
            "  height = %d,\n"
            "  transform = %s,\n"
            "]",
            m_resolution.x(),
            m_resolution.y(),
            indent(m_transform));
    }

private:
    /// @brief Returns the distance between the lens system and the film plane.
    float lensFrontZ() const {
        return m_elementInterfaces.front()->getThickness();
    }

    /// @brief Returns distance between lens system and the aperture plane.
    float lensRearZ() const {
        float zSum = 0;
        for (auto &element : m_elementInterfaces)
            zSum += element->getThickness();
        return zSum;
    }

    /// @brief Returns the radius of the front lens.
    float frontLensRadius() const {
        return m_elementInterfaces.front()->getApertureRadius();
    }

    /**
     * @brief Given a ray starting from the film side of the lens system,
     * computes intersections with each element in turn, terminating the ray
     * and returning false if its path is blocked along the way through the lens
     * system.
     */
    bool traceLensesFromFilm(const Ray &cameraRay, Ray &outRay) const {

        float elementZ = 0; // Start at the film plane
        Ray lensRay    = cameraRay;
        float etaT = 1; 
        for (size_t i = 0; i < m_elementInterfaces.size(); ++i) {

            const auto &element = m_elementInterfaces[i];
            elementZ += element->getThickness();

            // Compute intersection of ray with lens element
            float t;
            Vector n;
            const bool isStop =
                (element->getCurvatureRadius() == 0); // Check for aperture stop
            if (isStop) {
                t = (elementZ - lensRay.origin.z()) / lensRay.direction.z();
            } else {
                if (!element->intersect(lensRay, elementZ, t, n))
                    return false;
            }

            // Test intersection point against element aperture
            const Point pHit = lensRay(t);
            const float r2   = sqr(pHit.x()) + sqr(pHit.y());
            if (r2 > sqr(element->getApertureRadius()))
                return false;
            lensRay.origin = pHit;

            // Update ray path for element interface interaction
            if (!isStop) {
                const float etaI = element->getEta();
                // Calculate the direction of the refracted ray
                const Vector lensRayDir = lensRay.direction.normalized();
                const Vector refractedDir =
                    refract(-lensRayDir, n, etaI / etaT);
                etaT = etaI;
                if (refractedDir.isZero())
                    return false;
                lensRay.direction = refractedDir;
            }
        }
        outRay = lensRay;
        return true;
    }

    /**
     * @brief Given a ray starting from the scene facing side of the lens
     * system, computes intersections with each element in turn, terminating the
     * ray and returning false if its path is blocked along the way through the
     * lens system.
     */
    bool traceLensesFromScene(const Ray &cameraRay, Ray &outRay) const {

        float elementZ =
            lensRearZ(); // Start from the scene facing side of the lens system

        // Trace ray through lens system
        Ray lensRay(cameraRay);
        float etaI, etaT;
        for (int i = m_elementInterfaces.size() - 1; i >= 0; --i) {

            const auto &element = m_elementInterfaces[i];

            // Compute intersection of ray with lens element
            float t;
            Vector n;
            const bool isStop =
                (element->getCurvatureRadius() == 0); // Check for aperture stop
            if (isStop) {
                t = (elementZ - lensRay.origin.z()) / lensRay.direction.z();
            } else {
                if (!element->intersect(lensRay, elementZ, t, n))
                    return false;
            }

            // Test intersection point against element aperture
            const Point pHit = lensRay(t);
            const float r2   = sqr(pHit.x()) + sqr(pHit.y());
            if (r2 > sqr(element->getApertureRadius()))
                return false;
            lensRay.origin = pHit;

            // Update ray path for element interface interaction
            if (!isStop) {
                etaT = m_elementInterfaces[i]->getEta();
                if (i >= 1) {
                    etaI = m_elementInterfaces[i - 1]->getEta();
                } else {
                    etaI = 1;
                }
                // Calculate the direction of the refracted ray
                const Vector lensRayDir = lensRay.direction.normalized();
                const Vector refractedDir =
                    refract(-lensRayDir, n, etaI / etaT);
                if (refractedDir.isZero())
                    return false;
                lensRay.direction = refractedDir;
            }
            elementZ -= element->getThickness();
        }

        // Transform from lens system space back to camera space
        outRay = lensRay;
        return true;
    }

    /// @brief Compute the cardinal points of one side of the lens system.
    void computeCardinalPoints(const Ray &inRay, const Ray &outRay, float &pz,
                               float &fz) const {
        // Find focal point
        float tf = -outRay.origin.x() / outRay.direction.x();
        fz       = outRay(tf).z();
        // Find principal point
        float tp =
            (inRay.origin.x() - outRay.origin.x()) / outRay.direction.x();
        pz = outRay(tp).z();
    }

    /// @brief Compute thick lens approximation for the lens system.
    void computeThickLensApproximation(float pz[2], float fz[2]) const {
        // Find height x from optical axis for parallel rays
        float x = .001f * m_sensorSize.length();
        // Compute cardinal points for scene side of lens system
        Ray rScene = Ray(Point(x, 0, lensRearZ() + 1), Vector(0, 0, -1));
        Ray rFilm;
        traceLensesFromScene(rScene, rFilm);
        computeCardinalPoints(rScene, rFilm, pz[0], fz[0]);
        // Compute cardinal points for film side of lens system
        rFilm = Ray(Point(x, 0, lensFrontZ() - 1), Vector(0, 0, 1));
        traceLensesFromFilm(rFilm, rScene);
        computeCardinalPoints(rFilm, rScene, pz[1], fz[1]);
    }

    /**
     * @brief Compute the distance from the rear element to the film plane
     * that will cause the lens system to focus at a given distance.
     */
    float focusThickLens() {
        float pz[2], fz[2];
        computeThickLensApproximation(pz, fz);
        logger(EInfo, "pz1 = %f", pz[0]);
        logger(EInfo, "pz2 = %f", pz[1]);
        // Compute translation of lens, delta, to focus at m_focusDistance
        float f = fz[1] - pz[1];
        logger(EInfo, "f = %f", f);
        const float z = m_focusDistance;

        const float b     = z - pz[1] - pz[0];
        const float c     = z * pz[0] - pz[0] * pz[1] - f * (z - pz[0] + pz[1]);
        const float delta = -0.5f * (-b + sqrt(sqr(b) + 4 * c));

        logger(EInfo, "delta = %f", delta);
        return m_elementInterfaces.front()->getThickness() + delta;
    }

    /// @brief Draw a single lens interface of the lens system.
    std::string drawInterface(float *zPos, int idx) {
        ref<Lens> interface = m_elementInterfaces[idx];
        const float r       = 1000 * interface->getCurvatureRadius();
        const float h       = 1000 * interface->getApertureRadius();
        const float d       = 1000 * interface->getThickness();

        *zPos += d;

        std::string interfaceDrawn;
        // Handle aperture stop.
        if (r == 0) {
            return drawAperture(*zPos, h);
        }

        float zIts;
        bool sf; // Used so that arc is drawn correctly.

        interfaceZIts(interface, *zPos, &zIts, &sf);

        // Lens surface is axis symmetric and ends at a point (h,zIts).
        interfaceDrawn = tfm::format(
            "<path d=\"M %f %f A %f %f 0 0 %d %f %f \" \n"
            "      fill=\"none\" \n"
            "      stroke=\"blue\" \n"
            "      stroke-width=\"0.5\"/>\n",
            zIts,
            h,
            r,
            r,
            sf,
            zIts,
            -h);

        // Two interfaces form a lens.
        if (idx != 0) {
            ref<Lens> interfacePrev = m_elementInterfaces[idx - 1];
            const float hPrev       = 1000 * interfacePrev->getApertureRadius();
            const float etaPrev     = interfacePrev->getEta();
            if (etaPrev != 1) {
                float zItsPrev;
                interfaceZIts(interfacePrev, *zPos - d, &zItsPrev, nullptr);
                interfaceDrawn += drawLensClosing(h, hPrev, zIts, zItsPrev);
            }
        }

        return interfaceDrawn;
    }

    /// @brief Convenience function to find z-coord where lens surface ends.
    void interfaceZIts(ref<Lens> interface, float z, float *zIts, bool *sf) {
        const float r = 1000 * interface->getCurvatureRadius();
        const float h = 1000 * interface->getApertureRadius();

        const float zCenter = z + r;

        // Find intersections of a line of height h with the interface.
        const float b = -2 * zCenter;
        const float c = sqr(h) + sqr(zCenter) - sqr(r);

        float t0, t1;
        solveQuadratic(b, c, t0, t1);
        // Case 1: Convex interface
        if (r > 0) {
            *zIts = t0;
            if (sf != nullptr)
                *sf = 1;
        }
        // Case 2: Concave interface
        else {
            *zIts = t1;
            if (sf != nullptr)
                *sf = 0;
        }
    }

    std::string drawOpticalAxis() {
        return "<line x1=\"0\" y1=\"0\" x2=\"400\" y2=\"0\""
               " style=\"stroke:black;stroke-width:0.5\" />"s;
    }

    std::string drawLine(float x1, float y1, float x2, float y2,
                         std::string_view color = "blue") {
        return tfm::format(
            "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" "
            " stroke=\"%s\" "
            " stroke-width=\"0.5\"/>\n",
            x1,
            y1,
            x2,
            y2,
            color);
    }

    std::string drawLensClosing(float h, float hPrev, float z, float zPrev) {
        const float hMax = max(hPrev, h);
        const float hMin = min(hPrev, h);

        std::string lineSegments;
        lineSegments += drawLine(zPrev, hMax, z, hMax);
        lineSegments += drawLine(zPrev, -hMax, z, -hMax);
        if (hMax == hPrev) {
            lineSegments += drawLine(z, hMin, z, hMax);
            lineSegments += drawLine(z, -hMin, z, -hMax);
        } else {
            lineSegments += drawLine(zPrev, hMin, zPrev, hMax);
            lineSegments += drawLine(zPrev, -hMin, zPrev, -hMax);
        }
        return lineSegments;
    }

    std::string drawAperture(float z, float h) {
        std::string apertureDrawn;
        apertureDrawn += drawLine(z, h, z, 4 * h, "gray");
        apertureDrawn += drawLine(z, -h, z, -4 * h, "gray");
        return apertureDrawn;
    }

    std::string createHeader() {
        return "<svg width=\"1000\" height=\"1000\" "
               " viewBox=\"0 -125 500 250\" "
               "xmlns=\"http://www.w3.org/2000/svg\">\n";
    }
};
} // namespace lightwave

REGISTER_CAMERA(Realistic, "realistic")