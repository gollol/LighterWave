#include <lightwave.hpp>
 
namespace lightwave {

#define GEOINFO_VARIABLE_MAPPING {     \
    {"position",        GeometryInfoVariables::POSITION},          \
    {"normal",          GeometryInfoVariables::NORMAL},            \
    {"tangent",         GeometryInfoVariables::TANGENT},           \
    {"trueNormal",      GeometryInfoVariables::TRUE_NORMAL},       \
    {"incoming",        GeometryInfoVariables::INCOMING},          \
    {"parametric",      GeometryInfoVariables::PARAMETRIC},        \
    {"backfacing",      GeometryInfoVariables::BACKFACING},        \
    {"pointiness",      GeometryInfoVariables::POINTINESS},        \
    {"randomPerIsland", GeometryInfoVariables::RANDOM_PER_ISLAND}, \
}

class GeometryInfo : public Texture {

protected:
    enum GeometryInfoVariables {
        POSITION, NORMAL, TANGENT, TRUE_NORMAL, INCOMING, PARAMETRIC, BACKFACING, POINTINESS, RANDOM_PER_ISLAND 
    };

    GeometryInfoVariables m_variable;
    
public:

    GeometryInfo(const Properties &properties) {
        m_variable = properties.getEnum<GeometryInfoVariables>("variable", GEOINFO_VARIABLE_MAPPING);
    }

    Color evaluate(const Point2 &uv, const Context &cont) const override {
        const Intersection its = *cont.its;

        switch(m_variable){
            case POSITION:
                return Color(Vector(its.position));
            case NORMAL:
                return Color(its.shadingNormal);
            case TANGENT:
                return Color(its.tangent);
            case TRUE_NORMAL:
                return Color(its.geometryNormal);
            case INCOMING:
                return Color(cont.negRayDirection);
            case PARAMETRIC:
                // Blender uses the uv coordinates as if the mesh was subdivided -> bunch of uv triangles
                // Not sure how to achieve that - single triangles match.
                return Color(its.uv.x(), its.uv.y(), 0);
            case BACKFACING:
                return Color(cont.backFacing);
            default:
                assert_condition(false, logger(EError, "GeometryInfo doesn't yet support Color type %s",
                    enumToString<GeometryInfoVariables>(m_variable, GEOINFO_VARIABLE_MAPPING)););
        }

        return Color(0);
    }

    float scalar(const Point2 &uv, const Context &cont) const override {

        switch(m_variable){
            case BACKFACING:
                return cont.backFacing;
            default:
                assert_condition(false, logger(EError, "GeometryInfo doesn't yet support scalar type %s",
                    enumToString<GeometryInfoVariables>(m_variable, GEOINFO_VARIABLE_MAPPING)););
        }
        return 0;
    }

    std::string toString() const override {
        return tfm::format("GeometryInfo[\n"
                            " variable = %s,\n"
                           "]", indent(m_variable));
    }
};

} // namespace lightwave

REGISTER_TEXTURE(GeometryInfo, "geometry_info")
