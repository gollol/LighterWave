#include <lightwave.hpp>

namespace lightwave {

const size_t SIZE_PIXEL_WARP = 2048;  //Resize

Point2 primaryToUvSpace(const Point2 &primary) {
    return { primary.x(), safe_acos(2 * primary.y() - 1) * InvPi };
}

Point2 uvToPrimarySpace(const Point2 &uv) {
    return { uv.x(), cos(uv.y() * Pi) / 2 + 0.5f };
}

class Warp2DEnvMap final : public BackgroundLight {
    /// @brief The texture to use as background
    ref<Texture> m_texture;
    /// @brief An optional transform from local-to-world space
    ref<Transform> m_transform;
    /// @brief The calculated prefix sum of m_texure in x direction for every y
    std::vector<std::vector<float>> x_prefix_sum;
    /// @brief The calculated prefix sum of m_texure in y direction for the sum in the x direction, scaled by the cosine and the SIZE_PIXEL_WARP
    std::vector<float> y_prefix_sum;
    /// @brief This enables MIS compensation
    bool m_mis_compensation;
    /// @brief This value controlls the stregth of MIS compensation
    float m_mis_compensation_strength;

    //Samples the prefix vector returns the index, the warped position [0.0, 1.0] and the pdf
    std::pair<size_t ,float> sample_prefix_vector(const std::vector<float> &prefix_list, const float sample) const {
        float size_prefix = *(prefix_list.end()-1);
        float scaled = sample * size_prefix;
        auto pos = std::upper_bound(prefix_list.begin(), prefix_list.end(), scaled);
        size_t index = pos - prefix_list.begin();
        
        float end = *pos;
        float start;
        if (index == 0){
            start = 0.0;
        } else {
            start = *(pos-1);
        }

        float pixel_size = 1.0f / SIZE_PIXEL_WARP;
        float size_index = end - start;
        float relativ_pos = (scaled - start) / size_index;
        float warped_pos = pixel_size * relativ_pos + pixel_size * float(index); //TODO check if this is correct

        return std::make_pair(index ,warped_pos);
    }

    Point2 sample_point(Point2 sample) const {
        auto y_sample = sample_prefix_vector(y_prefix_sum, sample.y());
        auto x_sample = sample_prefix_vector(x_prefix_sum[std::get<0>(y_sample)], sample.x());
        Point2 warped_pos = Point2(std::get<1>(x_sample), std::get<1>(y_sample));
        return primaryToUvSpace(warped_pos);
    }

    inline float sample_pdf_pos(Point2 pos) const {
        pos = uvToPrimarySpace(pos);
        size_t x_pos = (size_t) min(pos.x() * float(SIZE_PIXEL_WARP), SIZE_PIXEL_WARP - 1); //this works because the implicit conversion from float to int is floor the min is just for the case pos.x == 1.0 where this might get unstable
        size_t y_pos = (size_t) min(pos.y() * float(SIZE_PIXEL_WARP), SIZE_PIXEL_WARP - 1);
        return sample_pdf_index(x_pos, y_pos);
    }

    inline float sample_pdf(float prefix_max_size, float size_index) const{
        if (prefix_max_size == 0.0f){
            return 0.0f;
        } else {
            return size_index * (float(SIZE_PIXEL_WARP) / (prefix_max_size));
        }
    }

    inline float sample_pdf_index(size_t x, size_t y) const {
        float x_prefix_size = x_prefix_sum[y][x_prefix_sum.size()-1];
        float x_start = x == 0 ? 0.0f : x_prefix_sum[y][x-1];
        float x_end = x_prefix_sum[y][x];
        float x_size = x_end - x_start;

        float y_prefix_size = y_prefix_sum[y_prefix_sum.size()-1];
        float y_start = y == 0 ? 0.0f : y_prefix_sum[y-1];
        float y_end = y_prefix_sum[y];
        float y_size = y_end - y_start;
        
        return sample_pdf(x_prefix_size, x_size) * sample_pdf(y_prefix_size, y_size);
    }

public:
    Warp2DEnvMap(const Properties &properties)
    : BackgroundLight(properties) {
        const Context cont = Context();
        m_texture   = properties.getChild<Texture>();
        m_transform = properties.getOptionalChild<Transform>();
        m_mis_compensation_strength = clamp(properties.get<float>("mis_compensation_strength"), 0.0f, 1.0f);
        m_mis_compensation = m_mis_compensation_strength > 0.001 ? true : false;

        //Calculate base image brightness and average brightness
        float one_pixel_size = 1.f / SIZE_PIXEL_WARP;
        float half_pixel_size = 0.5f * one_pixel_size; 
        std::vector<std::vector<float>> weights = std::vector<std::vector<float>>(SIZE_PIXEL_WARP, std::vector<float>(SIZE_PIXEL_WARP, 0.0f));
        float sum = 0.0f;
        for (size_t y = 0; y < SIZE_PIXEL_WARP; y++){
            float temp_sum = 0.0f; //We use this to avoid numeric problems 
            for (size_t x = 0; x < SIZE_PIXEL_WARP; x++){
                float x_pos = float(x) * one_pixel_size + half_pixel_size;
                float y_pos = float(y) * one_pixel_size + half_pixel_size;
                float mean = m_texture->evaluate(primaryToUvSpace(Point2(x_pos, y_pos)), cont).mean();
                weights[y][x] = mean;
                temp_sum += mean;
            }
            sum += temp_sum / SIZE_PIXEL_WARP;
        }
        float avr = sum / SIZE_PIXEL_WARP;

        //MIS compensation
        if (m_mis_compensation){
            float mis_constant = 2.0f * m_mis_compensation_strength * avr;
            for (size_t y = 0; y < SIZE_PIXEL_WARP; y++){
                for (size_t x = 0; x < SIZE_PIXEL_WARP; x++){
                    weights[y][x] = max(weights[y][x] - mis_constant, 0.0f);
                }
            }
        }

        //Precalculate X-Prefix Sum
        x_prefix_sum = std::vector<std::vector<float>>(SIZE_PIXEL_WARP, std::vector<float>(SIZE_PIXEL_WARP, 0.0));
        for (size_t y = 0; y < SIZE_PIXEL_WARP; y++){
            float prefix_sum = 0.0;
            for (size_t x = 0; x < SIZE_PIXEL_WARP; x++){
                float x_weight = weights[y][x];
                prefix_sum += x_weight;
                x_prefix_sum[y][x] = prefix_sum;
            }
        }

        //Precalculate Y-Prefix Sum
        y_prefix_sum = std::vector<float>(SIZE_PIXEL_WARP, 0.0);
        float prefix_sum = 0.0;
        for (size_t y = 0; y < SIZE_PIXEL_WARP; y++){
            float x_data = x_prefix_sum[y][SIZE_PIXEL_WARP-1];
            float y_pos = float(y) * one_pixel_size + half_pixel_size;
            float y_pos_angle = (y_pos - 0.5f) * Pi; 
            float t_x_data = x_data * abs(cos(y_pos_angle)) / float (SIZE_PIXEL_WARP); //We divide by the SIZE_PIXEL_WARP to avoid nummeric instabilty
            float y_weight = t_x_data;
            prefix_sum += y_weight;
            y_prefix_sum[y] = prefix_sum;
        }
    }

    DirectLightSample sampleDirect(const Point &origin,
                                   Sampler &rng) const override {
        PROFILE("Warp Direct")
        const Context cont = Context();
        const Point2 warped  = sample_point(rng.next2D());
        const float cosTheta = cos(warped.y() * Pi);
        const float sinTheta = safe_sqrt(1 - sqr(cosTheta));

        Vector direction = {
            sinTheta * -cos(2 * Pi * warped.x()),
            cosTheta,
            sinTheta * sin(2 * Pi * warped.x()),
        };

        if (m_transform) {
            SimpleTransform T = m_transform->interpolate(0);
            direction         = T.apply(direction).normalized();
        }

        const auto E    = m_texture->evaluate(warped, cont);
        const float pdf = sample_pdf_pos(warped) * Inv4Pi;
        
        return {
            .wi = direction,
            .weight = E / pdf,
            .pdf = pdf,
            .distance = Infinity,
        };
    }

    EmissionEval evaluate(const Vector &direction) const override {
        PROFILE("Warp Eval")
        const Context cont = Context();
        Vector d = direction;
        if (m_transform) {
            SimpleTransform T = m_transform->interpolate(0);
            d                 = T.inverse(d).normalized();
        }
        const Point2 warped((std::atan2(-d.z(), d.x()) + Pi) * Inv2Pi,
                            safe_acos(d.y()) * InvPi);
        
        assert_finite(m_texture->evaluate(warped, cont),
                { logger(EWarn, "warped: %s from %s", warped, d); });
        return {
            .value = m_texture->evaluate(warped, cont),
            .pdf = sample_pdf_pos(warped) * Inv4Pi,
        };
    }

    std::string toString() const override {
        return tfm::format(
            "Warp2DEnvMap[\n"
            "   texture = %s,\n"
            "   transform = %s,\n"
            "   m_mis_compensation = %s,\n"
            "   m_mis_compensation_strength = %s"
            "]",
            indent(m_texture),
            indent(m_transform),
            indent(m_mis_compensation),
            indent(m_mis_compensation_strength)
            );
    }
};

} // namespace lightwave

REGISTER_LIGHT(Warp2DEnvMap, "warp2devnmap")
