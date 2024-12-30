#include <lightwave.hpp>

namespace lightwave {

#define OINFO_VARIABLE_MAPPING {     \
    {"location",      ObjectInfoVariables::LOCATION},        \
    {"color",         ObjectInfoVariables::COLOR},           \
    {"alpha",         ObjectInfoVariables::ALPHA},           \
    {"objectIndex",   ObjectInfoVariables::OBJECT_INDEX},    \
    {"materialIndex", ObjectInfoVariables::MATERIAL_INDEX},  \
    {"random",        ObjectInfoVariables::RANDOM},          \
}

class ObjectInfo : public Texture {

public:

    enum ObjectInfoVariables {
        LOCATION, COLOR, ALPHA, OBJECT_INDEX, MATERIAL_INDEX, RANDOM
    };

    ObjectInfoVariables m_variable;

    ObjectInfo(const Properties &properties) {
        m_variable = properties.getEnum<ObjectInfoVariables>("variable", OINFO_VARIABLE_MAPPING);
    }

    Color evaluate(const Point2 &uv, const Context &cont) const override {
        switch (m_variable){
            case RANDOM:{
                int index = cont.its->instance->index();
                return Color(hash::hash_uint2_to_float(index, 0));
            }
            default:
                assert_condition(false, logger(EError, "ObjectInfo doesn't yet support color type %s",
                    enumToString<ObjectInfoVariables>(m_variable, OINFO_VARIABLE_MAPPING)););
        }

        return Color(0);
    }

    float scalar(const Point2 &uv, const Context &cont) const override {
        int index = cont.its->instance->index();

        switch (m_variable){
            case OBJECT_INDEX:
                return 1.f*index;
            case RANDOM:
                return hash::hash_uint2_to_float(index, 0);
            default:
                assert_condition(false, logger(EError, "ObjectInfo doesn't yet support scalar type %s",
                    enumToString<ObjectInfoVariables>(m_variable, OINFO_VARIABLE_MAPPING)););
        }
        return 0;
    }

    std::string toString() const override {
        return tfm::format(
                        "ObjectInfo[\n"
                        "  variable = %s,\n"
                        "]", indent(m_variable));
    }
};

} // namespace lightwave

REGISTER_TEXTURE(ObjectInfo, "object_info")
