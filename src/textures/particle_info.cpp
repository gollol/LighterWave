#include <lightwave.hpp>
 
namespace lightwave {

#define PINFO_VARIABLE_MAPPING {     \
    {"index",           ParticleInfoVariables::INDEX},             \
    {"random",          ParticleInfoVariables::RANDOM},            \
    {"age",             ParticleInfoVariables::AGE},               \
    {"size",            ParticleInfoVariables::SIZE},              \
    {"lifetime",        ParticleInfoVariables::LIFETIME},          \
    {"location",        ParticleInfoVariables::LOCATION},          \
    {"velocity",        ParticleInfoVariables::VELOCITY},          \
    {"angularVelocity", ParticleInfoVariables::ANGULAR_VELOCITY},  \
}

class ParticleInfo : public Texture {

protected:
    enum ParticleInfoVariables {
        INDEX, RANDOM, AGE, LIFETIME, LOCATION, SIZE, VELOCITY, ANGULAR_VELOCITY
    };

    ParticleInfoVariables m_variable;
    
public:

    ParticleInfo(const Properties &properties) {
        m_variable = properties.getEnum<ParticleInfoVariables>("variable", PINFO_VARIABLE_MAPPING);      
    }

    // We don't evaluate Color from particle Infos!
    Color evaluate(const Point2 &uv, const Context &cont) const override {
        return Color(scalar(uv, cont));
    }

    float scalar(const Point2 &uv, const Context &cont) const override {
        int index = max(cont.its->instance->index()-1, -1);   // Particles now start at 1 instead of 0, non-particles used to be -1.
        
        switch(m_variable){
            case INDEX:
                return 1.f*index;
            case RANDOM:
                return hash::hash_uint2_to_float(index, 0);
            default:
                assert_condition(false, logger(EError, "ParticleInfo doesn't yet support type %s",
                    enumToString<ParticleInfoVariables>(m_variable, PINFO_VARIABLE_MAPPING)););
        }
        return 0;
    }

    std::string toString() const override {
        return tfm::format("ParticleInfo[\n"
                            " variable = %s,\n"
                           "]", indent(m_variable));
    }
};

} // namespace lightwave

REGISTER_TEXTURE(ParticleInfo, "particle_info")
