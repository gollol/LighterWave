#include <lightwave.hpp>

namespace lightwave {

/*
 * These defines are used, so we can convert the enums back to strings
*/

#define DIMENSIONS_MAPPING {        \
        { "1D", Dimensions::D1 },   \
        { "2D", Dimensions::D2 },   \
        { "3D", Dimensions::D3 },   \
        { "4D", Dimensions::D4 },   \
    }

#define NOISE_TYPE_MAPPING {                                        \
        { "FBM", NoiseType::FBM },                                  \
        { "MULTIFRACTAL", NoiseType::MULTIFRACTAL },                \
        { "HYBRID_MULTIFRACTAL", NoiseType::HYBRID_MULTIFRACTAL },  \
        { "RIDGED_MULTIFRACTAL", NoiseType::RIDGED_MULTIFRACTAL },  \
        { "HETERO_TERRAIN", NoiseType::HETERO_TERRAIN },            \
    }

class NoiseTexture : public Texture {

    // Decribes the available dimensionalities for the input to the noise texture
    /// @note Currently, only 2D is supported, because we only have UV coordinates
    enum class Dimensions {
        D1,
        D2,
        D3,
        D4,
    };

    // Type of Noise texture, with different ways to combine octaves
    /// @note Currently, only fBM is supported  
    enum class NoiseType {
        FBM,
        MULTIFRACTAL,
        HYBRID_MULTIFRACTAL,
        RIDGED_MULTIFRACTAL,
        HETERO_TERRAIN,
    };

    ref<Texture> m_scale;       // default: 5
    ref<Texture> m_detail;      // default: 2
    ref<Texture> m_roughness;   // default: 0.5
    ref<Texture> m_lacunarity;  // default: 2
    ref<Texture> m_distortion;  // default: 0

    Dimensions m_dimensions;
    NoiseType m_noiseType;
    bool m_normalize;
    bool m_colorful;

    Vector2 randomVector2Offset(float seed) const {
        using namespace lightwave::hash;

        return Vector2(
            100.0f + hash_Point2_to_float(Point2(seed, 0.0f)) * 100.0f,
            100.0f + hash_Point2_to_float(Point2(seed, 1.0f)) * 100.0f
        );
    }

    /// @note See https://projects.blender.org/blender/blender/src/commit/4651f8a08faa578ff1a488e68fe8ee659cd31424/intern/cycles/kernel/svm/noisetex.h#L118
    Color noiseTexture2D(Point2 p, const float detail, const float roughness,
                         const float lacunarity, const float distortion, const bool normalize, 
                         const bool colorful) const {
        if (distortion != 0.0f) {
            // We have distortion
            p += Vector2(
                snoise_2d(p + randomVector2Offset(0.0f)) * distortion,
                snoise_2d(p + randomVector2Offset(1.0f)) * distortion
            );
        }

        float value = noise_fbm(p, detail, roughness, lacunarity, normalize);

        if (colorful) {
            return Color(
                value,
                noise_fbm(p + randomVector2Offset(2.0f), detail, roughness, lacunarity, normalize),
                noise_fbm(p + randomVector2Offset(3.0f), detail, roughness, lacunarity, normalize)
            );
        } else {
            return Color(value);
        }
    }

public:
    NoiseTexture(const Properties &properties) {
        this->m_scale  = properties.get<Texture>("scale");
        this->m_detail  = properties.get<Texture>("detail");
        this->m_roughness  = properties.get<Texture>("roughness");
        this->m_lacunarity  = properties.get<Texture>("lacunarity");
        this->m_distortion  = properties.get<Texture>("distortion");

        this->m_normalize  = properties.get<bool>("normalize", true);
        this->m_colorful  = properties.get<bool>("colorful", false);

        this->m_dimensions = properties.getEnum<Dimensions>("dimensions", Dimensions::D3, DIMENSIONS_MAPPING);

        assert_condition(this->m_dimensions == Dimensions::D2, {
            logger(EError, "The specified dimensionality is currently not supported.");
        });

        this->m_noiseType = properties.getEnum<NoiseType>("noisetype", NoiseType::FBM, NOISE_TYPE_MAPPING);

        assert_condition(this->m_noiseType == NoiseType::FBM, {
            logger(EError, "The specified Noise Texture Type is currently not supported. %s", this);
        });
    }

    Color evaluate(const Point2 &uv, const Context &cont) const override {
        return noiseTexture2D(
            Point2(uv.x() * this->m_scale->scalar(uv, cont), uv.y() * this->m_scale->scalar(uv, cont)),
            this->m_detail->scalar(uv, cont),
            this->m_roughness->scalar(uv, cont),
            this->m_lacunarity->scalar(uv, cont),
            this->m_distortion->scalar(uv, cont),
            this->m_normalize,
            this->m_colorful
        );
    }

    std::string toString() const override {
        return tfm::format("NoiseTexture[\n"
                           "  scale = %f\n"
                           "  detail = %f\n"
                           "  roughness = %f\n"
                           "  lacunarity = %f\n"
                           "  distortion = %f\n"
                           "  normalize = %s\n"
                           "  colorful = %s\n"
                           "  dimensions = %s\n"
                           "  type = %s\n"
                           "]", 
                           indent(this->m_scale),
                           indent(this->m_detail),
                           indent(this->m_roughness),
                           indent(this->m_lacunarity),
                           indent(this->m_distortion),
                           this->m_normalize,
                           this->m_colorful,
                           enumToString<Dimensions>(this->m_dimensions, DIMENSIONS_MAPPING),
                           enumToString<NoiseType>(this->m_noiseType, NOISE_TYPE_MAPPING)
                           );
    }
};

} // namespace lightwave

REGISTER_TEXTURE(NoiseTexture, "noise")
