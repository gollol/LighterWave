#include <lightwave.hpp>

namespace lightwave {

class AgX : public Postprocess {
//https://iolite-engine.com/blog_posts/minimal_agx_implementation
public:
    AgX(const Properties &properties) : Postprocess(properties) 
    {
    }

    Color agxDefaultContrastApprox(Color x) {
        //Color x2 = x * x;
        //Color x4 = x2 * x2;
    //
        //return  + 15.5     * x4 * x2
        //        - 40.14    * x4 * x
        //        + 31.96    * x4
        //        - 6.868    * x2 * x
        //        + 0.4298   * x2
        //        + 0.1191   * x
        //        - 0.00232  * Color(1);

        Color x2 = x * x;
        Color x4 = x2 * x2;
        Color x6 = x4 * x2;

        return - 17.86     * x6 * x
               + 78.01     * x6
               - 126.7     * x4 * x
               + 92.06     * x4
               - 28.72     * x2 * x
               + 4.361     * x2
               - 0.1718    * x
               + 0.002857  * Color(1);
    }

    Color agx(Color val) {
        const Matrix3x3 agx_mat = Matrix3x3({
          0.842479062253094, 0.0423282422610123, 0.0423756549057051,
          0.0784335999999992,  0.878468636469772,  0.0784336,
          0.0792237451477643, 0.0791661274605434, 0.879142973793104}).transpose();

        const float min_ev = -12.47393f;
        const float max_ev = 4.026069f;

        // Input transform (inset)
        val = Color(agx_mat * Vector(val.r(), val.g(), val.b()));

        // Log2 space encoding
        val.r() = clamp(log2(val.r()), min_ev, max_ev);
        val.g() = clamp(log2(val.g()), min_ev, max_ev);
        val.b() = clamp(log2(val.b()), min_ev, max_ev);

        val = (val - Color(min_ev)) / (max_ev - min_ev);

        // Apply sigmoid function approximation
        val = agxDefaultContrastApprox(val);

        return val;
    }

    Color agxEotf(Color val) {
        const Matrix3x3 agx_mat_inv = Matrix3x3({
            1.19687900512017, -0.0528968517574562, -0.0529716355144438,
            -0.0980208811401368, 1.15190312990417, -0.0980434501171241,
            -0.0990297440797205, -0.0989611768448433, 1.15107367264116}).transpose();

        // Inverse input transform (outset)
        val = Color(agx_mat_inv * Vector(val.r(), val.g(), val.b()));

        // sRGB IEC 61966-2-1 2.2 Exponent Reference EOTF Display
        // NOTE: We're linearizing the output here. Comment/adjust when
        // *not* using a sRGB render target
        val.r() = pow(val.r(), 2.2f);
        val.g() = pow(val.g(), 2.2f);
        val.b() = pow(val.b(), 2.2f);

        return val;
    }

    Color agxLook(Color val) {
        const Color lw = Color(0.2126, 0.7152, 0.0722);
        //float luma = dot(val, lw);
        float luma = val.r() * lw.r() + val.g() * lw.g() + val.b() + lw.b();

        // Default
        Color offset = Color(0.0);
        Color slope = Color(1.0);
        Color power = Color(1.0);
        float sat = 1.0;
          
        //#define AGX_LOOK 0

        #if AGX_LOOK == 1
          // Golden
          slope = Color(1.0, 0.9, 0.5);
          power = Color(0.8);
          sat = 0.8;
        #elif AGX_LOOK == 2
          // Punchy
          slope = Color(1.0);
          power = Color(1.35, 1.35, 1.35);
          sat = 1.4;
        #endif

          // ASC CDL
          val.r() = pow(val.r() * slope.r() + offset.r(), power.r());
          val.g() = pow(val.g() * slope.g() + offset.g(), power.g());
          val.b() = pow(val.b() * slope.b() + offset.b(), power.b());

          return Color(luma) + sat * (val - Color(luma));
    }

    /// @brief AgX nach https://iolite-engine.com/blog_posts/minimal_agx_implementation
    void execute() override 
    {
        //Output initialized to resolution of input
        m_output->initialize(m_input->resolution());
        
        for(int y = 0; y < m_output->resolution().y(); y++)
        for(int x = 0; x < m_output->resolution().x(); x++)
        {
            Point2i pixel(x,y);
            Color rgb = m_input->get(pixel);

            rgb = agx(rgb);
            //rgb = agxLook(rgb);
            rgb = agxEotf(rgb);

            m_output->get(pixel) = rgb;
        }

        m_output->save();        
    }

    std::string toString() const override {
        return tfm::format(
            "AgX[\n"
            "  input = %s,\n"
            "  output = %s,\n"
            "]",
            indent(m_input),
            indent(m_output)
        );
    }
};

} // namespace lightwave

REGISTER_POSTPROCESS(AgX, "agx");
