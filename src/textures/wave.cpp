#include <lightwave.hpp>

namespace lightwave {

#define WAVE_TYPE_MAPPING {             \
        { "BANDS", WaveType::BANDS },   \
        { "RINGS", WaveType::RINGS },   \
    }

#define BANDS_DIRECTION_MAPPING {                   \
        { "X", BandsDirection::X },                 \
        { "Y", BandsDirection::Y },                 \
        { "Z", BandsDirection::Z },                 \
        { "DIAGONAL", BandsDirection::DIAGONAL },   \
    }

#define WAVE_PROFILE_MAPPING {              \
        { "SIN", WaveProfile::SINE },       \
        { "SAW", WaveProfile::SAW },        \
        { "TRI", WaveProfile::TRIANGLE },   \
    }

class WaveTexture : public Texture {

    // Decribes the type of wave to be generated
    /// @note Currently, only BANDS is supported
    enum class WaveType {
        BANDS,
        RINGS,
    };

    // In which direction the wave bands should go
    /// @note Currently, Z is not supported, as we only have 2D UV textures  
    enum class BandsDirection {
        X,
        Y,
        Z,
        DIAGONAL,
    };

    // Describes what profile the waves should follow
    enum class WaveProfile {
        SINE,
        SAW,
        TRIANGLE,
    };

    ref<Texture> m_scale;           // default: 5
    ref<Texture> m_distortion;      // default: 0
    ref<Texture> m_detail;          // default: 2
    ref<Texture> m_detailScale;     // default: 1 
    ref<Texture> m_detailRoughness; // default: 0.5
    ref<Texture> m_phaseOffset;     // default: 0

    WaveType m_waveType;
    BandsDirection m_bandsDirection;
    WaveProfile m_waveProfile;

    /// @note See https://projects.blender.org/blender/blender/src/commit/5ce29bedf63054c0487913ae2b15f7617d95d1f7/intern/cycles/kernel/svm/wave.h
    float waveBandsTexture2D(Point2 p, const BandsDirection direction, const WaveProfile profile,
                             const float distortion, const float detail, const float detailScale,
                             const float detailRoughness, const float phaseOffset) const {
        //p = (p + Vector2(0.000001f)) * 0.999999f; // Used by Blender, but still not fully faithfull

        float n;

        // Add selected coordinate
        switch (direction) {
        case BandsDirection::X:
            n = p.x() * 20.0f;
            break;
        case BandsDirection::Y:
            n = p.y() * 20.0f;
            break;
        case BandsDirection::Z:
            // This case is not implemented as UV coordinates are only 2D
            // It is already caught by the constructor
            break;
        case BandsDirection::DIAGONAL:
            n = (p.x() + p.y()) * 10.0f;
            break;
        }

        // Add phase offset and distortion using FBM noise
        n += phaseOffset;

        if (distortion != 0.0f) {
            n += distortion * (noise_fbm(p * detailScale, detail, detailRoughness, 2.0f, true) * 2.0f - 1.0f);
        }

        switch (profile) {
        case WaveProfile::SINE:
            return 0.5f + 0.5f * sinf(n - Pi2);
        case WaveProfile::SAW:
            n /= 2 * Pi;
            return n - floorf(n);
        case WaveProfile::TRIANGLE:
            n /= 2 * Pi;
            return fabsf(n - floorf(n + 0.5f)) * 2.0f;
        }
    }

public:
    WaveTexture(const Properties &properties) {
        this->m_scale  = properties.get<Texture>("scale");
        this->m_distortion  = properties.get<Texture>("distortion");
        this->m_detail  = properties.get<Texture>("detail");
        this->m_detailScale  = properties.get<Texture>("detailscale");
        this->m_detailRoughness  = properties.get<Texture>("detailroughness");
        this->m_phaseOffset  = properties.get<Texture>("phaseoffset");

        this->m_waveType = properties.getEnum<WaveType>("wavetype", WaveType::BANDS, WAVE_TYPE_MAPPING);

        assert_condition(this->m_waveType == WaveType::BANDS, {
            logger(EError, "The specified wave type '%s' is currently not supported.", enumToString<WaveType>(this->m_waveType, WAVE_TYPE_MAPPING));
        });

        this->m_bandsDirection = properties.getEnum<BandsDirection>("direction", BandsDirection::X, BANDS_DIRECTION_MAPPING);

        assert_condition(this->m_bandsDirection != BandsDirection::Z, {
            logger(EError, "The bands direction Z is not supported, as lighwave only has 2D UV coordinates.");
        });

        this->m_waveProfile = properties.getEnum<WaveProfile>("waveprofile", WaveProfile::SINE, WAVE_PROFILE_MAPPING);
    }

    Color evaluate(const Point2 &uv, const Context &cont) const override {
        switch (this->m_waveType) {
        case WaveType::BANDS:
            return Color(waveBandsTexture2D(
                Point2(uv.x() * this->m_scale->scalar(uv, cont), uv.y() * this->m_scale->scalar(uv, cont)),
                this->m_bandsDirection,
                this->m_waveProfile,
                this->m_distortion->scalar(uv, cont),
                this->m_detail->scalar(uv, cont),
                this->m_detailScale->scalar(uv, cont),
                this->m_detailRoughness->scalar(uv, cont),
                this->m_phaseOffset->scalar(uv, cont)
            ));
        case WaveType::RINGS:
            // Currently not supported
            return Color(0);
        }
    }

    std::string toString() const override {
        return tfm::format("WaveTexture[\n"
                           "  waveType = %s\n"
                           "  bandsDirection = %s\n"
                           "  waveProfile = %s\n"
                           "  scale = %s\n"
                           "  distortion = %s\n"
                           "  detail = %s\n"
                           "  detailScale = %s\n"
                           "  detailRoughness = %s\n"
                           "  phaseOffset = %s\n"
                           "]", 
                           enumToString<WaveType>(this->m_waveType, WAVE_TYPE_MAPPING),
                           enumToString<BandsDirection>(this->m_bandsDirection, BANDS_DIRECTION_MAPPING),
                           enumToString<WaveProfile>(this->m_waveProfile, WAVE_PROFILE_MAPPING),
                           indent(this->m_scale),
                           indent(this->m_distortion),
                           indent(this->m_detail),
                           indent(this->m_detailScale),
                           indent(this->m_detailRoughness),
                           indent(this->m_phaseOffset)
                           );
    }
};

} // namespace lightwave

REGISTER_TEXTURE(WaveTexture, "wave")
