#include <lightwave.hpp>
 
namespace lightwave {

#define CINFO_VARIABLE_MAPPING {     \
    {"isStrand",      CurveInfoVariables::IS_STRAND},        \
    {"intercept",     CurveInfoVariables::INTERCEPT},        \
    {"length",        CurveInfoVariables::LENGTH},           \
    {"thickness",     CurveInfoVariables::THICKNESS},        \
    {"tangentNormal", CurveInfoVariables::TANGENT_NORMAL},   \
    {"random",        CurveInfoVariables::RANDOM},           \
}

// Because of non-uniform hair strand length, we don't support the length variable yet!
class CurvesInfo : public Texture {

protected:
    enum CurveInfoVariables {
        IS_STRAND, INTERCEPT, LENGTH, THICKNESS, TANGENT_NORMAL, RANDOM
    };

    CurveInfoVariables m_variable;
    
public:
    CurvesInfo(const Properties &properties) {
        m_variable = properties.getEnum<CurveInfoVariables>("variable", CINFO_VARIABLE_MAPPING);
    }

    Color evaluate(const Point2 &uv, const Context &cont) const override {
        if(m_variable == TANGENT_NORMAL)
            return Color(cont.tangentNormal);
        return Color(scalar(uv, cont));
    }

    float scalar(const Point2 &uv, const Context &cont) const override {
        int curveIndex = cont.curvesIndex;

        switch(m_variable){
            case IS_STRAND:
                return curveIndex == -1 ? 0 : 1;
            case INTERCEPT:
                return cont.hairIntercept;
            case THICKNESS:
                return cont.hairThickness;
            case RANDOM:{
                if(curveIndex == -1)
                    return 0;
                return hash::hash_uint2_to_float(curveIndex, 0);
            }
            default:
                assert_condition(false, logger(EError, "CurvesInfo doesn't yet support type %s",
                    enumToString<CurveInfoVariables>(m_variable, CINFO_VARIABLE_MAPPING)););
        }
        return 0;
    }

    std::string toString() const override {
        return tfm::format("CurvesInfo[\n"
                            " variable = %s,\n"
                           "]", indent(m_variable));
    }
};

} // namespace lightwave

REGISTER_TEXTURE(CurvesInfo, "curves_info")
