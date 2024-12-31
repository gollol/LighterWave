#include <lightwave.hpp>

namespace lightwave {

class Tonemapping : public Postprocess {
    float epsilson = 0.001f*100;

    Image r_m;

    Image gauss_image;

    Image gauss1;
    Image gauss2;

    float a;

enum Variable {
    Draco,
    ReinhardSimple,
    ReinhardLocal,
} m_operation;

public:
    Tonemapping(const Properties &properties) : Postprocess(properties) 
    {
        a = properties.get<float>("key");
        m_operation = properties.getEnum<Variable>("tmmode", {
            { "ReinhardSimple",  ReinhardSimple},
            { "ReinhardLocal",  ReinhardLocal},
            { "Draco",  Draco},
        });
    }

    /// @brief Draco TM
    void draco() 
    {
        m_output->initialize(m_input->resolution());

        float sum_log = 0.0f;
        float L_max = -Infinity;
        for(int c = 0; c < 3; c++)
        {
            for(int y = 0; y < m_output->resolution().y(); y++)
            for(int x = 0; x < m_output->resolution().x(); x++)
            {
                Point2i pixel(x,y);

                Color rgb = m_input->get(pixel);
                sum_log += log(rgb[c] + epsilson);  
                if(L_max < rgb[c])
                    L_max = rgb[c];
            }

            for(int y = 0; y < m_output->resolution().y(); y++)
            for(int x = 0; x < m_output->resolution().x(); x++)
            {
                Point2i pixel(x,y);
                Color rgb = m_input->get(pixel);

                //Drago Log TM
                float bias = log(0.8f) / log(0.5f);
                float L = log(1.0f + rgb[c]) / (log10(1.0f + L_max) * log(2.0f + 8.0f * pow(rgb[c] / L_max, bias)));

                float gamma = 2.2f;
                L = pow(L, 1.0f / gamma);

                m_output->get(pixel)[c] = L;
            }

            sum_log = 0.0f;
            L_max = -Infinity;
        }

        m_output->save();
    }

    /// @brief Reinhard Simple TM
    void reinhardSimple()
    {
        //init
        m_output->initialize(m_input->resolution());
        r_m.initialize(m_input->resolution());

        //user specified values
        // float a = 0.18f;
        // float phi = 10.0f;

        //calcute intermediate image r_m
        float sum_log = 0.0f;
        for(int c = 0; c < 3; c++)
        {
            //calcute perceived average luminance r_lav.
            for(int y = 0; y < m_output->resolution().y(); y++)
            for(int x = 0; x < m_output->resolution().x(); x++)
            {
                Point2i pixel(x,y);

                Color rgb = m_input->get(pixel);
                sum_log += log(rgb[c] + epsilson*1);  
            }
            float r_lav = exp(sum_log / (m_output->resolution().x() * m_output->resolution().y()) -  epsilson*1);

            //build r_m
            for(int y = 0; y < m_output->resolution().y(); y++)
            for(int x = 0; x < m_output->resolution().x(); x++)
            {
                Point2i pixel(x,y);
                Color rgb = m_input->get(pixel);

                //float r_m = a / r_lav * rgb[c];
                r_m.get(pixel)[c] = a / r_lav * rgb[c];

                m_output->get(pixel)[c] = r_m.get(pixel)[c] / (1.0f + r_m.get(pixel)[c]);
            }

            sum_log = 0.0f;
        }

        m_output->save();
    }

    std::vector<float> buildKernel(int size)
    {
        float c = 1.0f / size;
        const int kernelSize = 2 * size + 1;
        std::vector<float> kernel(kernelSize, 0.f);

        float kernel_sum = 0.0f;
        for(int i = 0; i < kernelSize; i++)
        {
            int x = i - size;
            //https://de.wikipedia.org/wiki/Gau%C3%9F-Filter
            kernel[i] = sqrt(c / Pi) * exp(-c*x*x);
            kernel_sum += kernel[i];
        }

        //make Kernel sum up to 1
        for(int i = 0; i < kernelSize; i++)
            kernel[i] /= kernel_sum;

        return kernel;
    }

    Color gauss_con(const std::vector<float> &kernel, const int size, Point2i pixel)
    {
        Color gauss_value = Color(0);
        const int kernelSize = static_cast<int>(kernel.size());

        const int sizePlus1 = size + 1;
        for(int y = pixel.y() - size; y <= pixel.y() + size; y++)
        for(int x = pixel.x() - size; x <= pixel.x() + size; x++)
        {
            int yy = y;
            if(y < 0)
                yy = 0;
            if(y >= r_m.resolution().y())
                yy = r_m.resolution().y() - 1;

            Color gauss_sum = Color(0.0f);

            for(int i = 0; i < kernelSize; i++)
            {
                int position = x + (i - size);

                if(position < 0) {
                    //position = -position;
                    position = 0;
                }

                if(position >= r_m.resolution().x())
                {
                    //position = img.resolution().x() - size;
                    position = r_m.resolution().x() - 1;
                }

                Color rgb = r_m.get(Point2i(position, yy));
                gauss_sum += kernel[i] * rgb;
            }

            Point2i pixel_gauss = Point2i(pixel.x() - (x - sizePlus1) - 1, pixel.y() - (yy - sizePlus1) - 1);
            gauss_image.get(pixel_gauss) = gauss_sum;
        }

        for(int y = 0; y < kernelSize; y++)
        {
            Color rgb = gauss_image.get(Point2i(sizePlus1, y));
            gauss_value += kernel[y] * rgb;
        }

        return gauss_value;
    }

    // This shit is never called???

    // void gauss_convolution(const int size, Image& img, Image& r_m)
    // {
    //     img.initialize(r_m.resolution());
        
    //     float c = 1.0f / size;

    //     const int kernelSize = 2 * size + 1;
    //     std::vector<float> kernel(kernelSize, 0.f);

    //     float kernel_sum = 0.0f;
    //     for(int i = 0; i < kernelSize; i++)
    //     {
    //         int x = i - size;
    //         //https://de.wikipedia.org/wiki/Gau%C3%9F-Filter
    //         kernel[i] = sqrt(c * InvPi) * exp(-c*x*x);
    //         kernel_sum += kernel[i];
    //     }

    //     for(int y = 0; y < img.resolution().y(); y++)
    //     for(int x = 0; x < img.resolution().x(); x++)
    //     {
    //         Color gauss_value = Color(0);
    //         for(int i = 0; i < kernelSize; i++)
    //         {
    //             int position = x + (i - size);

    //             if(position < 0) {
    //                 //position = -position;
    //                 position = 0;
    //             }

    //             if(position >= img.resolution().x()) {
    //                 //position = img.resolution().x() - size;
    //                 position = img.resolution().x() - 1;
    //             }

    //             Color rgb = r_m.get(Point2i(position, y));
    //             gauss_value += kernel[i] * rgb;
    //         }
    //         gauss_value /= kernel_sum;

    //         Point2i pixel(x,y);
    //         img.get(pixel) = gauss_value;
    //     }

    //     for(int y = 0; y < img.resolution().y(); y++)
    //     for(int x = 0; x < img.resolution().x(); x++)
    //     {
    //         Color gauss_value = Color(0);
    //         for(int i = 0; i < kernelSize; i++)
    //         {
    //             int position = y + (i - size);

    //             if(position < 0)
    //                 //position = -position;
    //                 position = 0;

    //             if(position >= img.resolution().y())
    //                 //position = img.resolution().x() - size;
    //                 position = img.resolution().y() - 1;

    //             Point2i pixel(x, position);
    //             Color rgb = img.get(pixel);

    //             gauss_value += kernel[i] * rgb;
    //         }
    //         gauss_value /= kernel_sum;

    //         Point2i pixel(x,y);
    //         img.get(pixel) = gauss_value;
    //     }
    // }

    void reinhardLocal()
    {
        //init
        m_output->initialize(m_input->resolution());
        r_m.initialize(m_input->resolution());

        //user specified values
        // float a = 0.18f;
        float phi = 10.0f;

        //calcute intermediate image r_m
        float sum_log = 0.0f;
        for(int c = 0; c < 3; c++)
        {
            //calcute perceived average luminance r_lav.
            for(int y = 0; y < m_output->resolution().y(); y++)
            for(int x = 0; x < m_output->resolution().x(); x++)
            {
                Color rgb = m_input->get(Point2i(x,y));
                sum_log += log(rgb[c] + epsilson*1);  
            }
            float r_lav = exp(sum_log / (m_output->resolution().x() * m_output->resolution().y()) -  epsilson*1);

            //build r_m
            for(int y = 0; y < m_output->resolution().y(); y++)
            for(int x = 0; x < m_output->resolution().x(); x++)
            {
                Point2i pixel(x,y);
                Color rgb = m_input->get(pixel);

                //float r_m = a / r_lav * rgb[c];
                r_m.get(pixel)[c] = a / r_lav * rgb[c];
            }

            sum_log = 0.0f;
        }
        
        //look for best gaussblurred image to adapt image locally (Dodging and Burning)
        for(int y = 0; y < m_output->resolution().y(); y++)
        for(int x = 0; x < m_output->resolution().x(); x++)
        {
            Point2i pixel(x,y);

            for(int c = 0; c < 3; c++)
            {
                int best_size = 0;
                float best_value = Infinity;
                for(int i = 1; i <= 3; i++)//TODO instead of 3 it can be user specified
                {
                    int size = i;
                    gauss_image.initialize(Point2i(2 * size + 1, 2 * size + 1));
                    std::vector<float> kernel = buildKernel(size);
                    Color g1 = gauss_con(kernel, size, pixel);

                    size++;
                    gauss_image.initialize(Point2i(2 * size + 1, 2 * size + 1));
                    kernel = buildKernel(size);
                    Color g2 = gauss_con(kernel, size, pixel);

                    float V = (g1 - g2)[c] / (pow(2.0f, phi) * (a / (i*i)) + g1[c]);
                    float VV = V * V;

                    if(abs(epsilson - VV) < best_value)
                    {
                        best_size = i;
                        best_value = abs(epsilson - VV);
                    }
                }

                gauss_image.initialize(Point2i(2 * best_size + 1, 2 * best_size + 1));
                std::vector<float> kernel = buildKernel(best_size);
                Color best_gauss = gauss_con(kernel, best_size, pixel);//You get a somekind of tape effect when using the original image and not r_m in gauss_con

                float T = r_m.get(pixel)[c] / (1.0f + best_gauss[c]);

                T = clamp(T, 0.0f, 1.0f);

                m_output->get(pixel)[c] = T;
            }
        }

        m_output->save();
    }

    /// @brief Tonemapping of image with local Reinhard Tonemapping
    void execute() override 
    {
        switch (m_operation)
        {
        case Draco:
            draco();
            break;
        case ReinhardSimple:
            reinhardSimple();
            break;
        case ReinhardLocal:
            reinhardLocal();
            break;
        default:
            draco();
            break;
        }
    }

    std::string toString() const override {
        return tfm::format(
            "Tonemapping[\n"
            "  input = %s,\n"
            "  output = %s,\n"
            "  a = %s, \n"
            "  TM type = %s, \n"
            "]",
            indent(m_input),
            indent(m_output),
            indent(a),
            indent(m_operation)
        );
    }
};

} // namespace lightwave

REGISTER_POSTPROCESS(Tonemapping, "tonemapping");

// Some color space transformations (HSV, YCbCr)

//    //https://gist.github.com/yohhoy/dafa5a47dade85d8b40625261af3776a with BT.601
//    inline Color rgb2yCbCr(Color rgb)
//    {
//        //float a = 0.299f;
//        //float b = 0.587f;
//        //float c = 0.114f;
//        //float d = 1.772f;
//        //float e = 1.402f;
////
//        //float Y = a * rgb.r() + b * rgb.g() + c * rgb.b();
////
//        //return Color(
//        //    Y, 
//        //    (rgb.b() - Y) / d, 
//        //    (rgb.r() - Y) / e
//        //);
//
//        return Color(
//                0.2627f         * rgb.r() + 0.678f      * rgb.g() + 0.0593f * rgb.b(),
//                -0.13963f       * rgb.r() - 0.360369f   * rgb.g() + 0.49999f * rgb.b(),
//                0.49999f        * rgb.r() - 0.459785f   * rgb.g() - 0.0402142f * rgb.b()
//        );
//    }
//
//    inline Color yCbCr2rgb(Color yCbCr)
//    {
//        //float a = 0.299f;
//        //float b = 0.587f;
//        //float c = 0.114f;
//        //float d = 1.772f;
//        //float e = 1.402f;
////
//        //return Color(
//        //    yCbCr.r() + e * yCbCr.b(), 
//        //    yCbCr.r() - (a * e / b) * yCbCr.g() - (c * d / b) * yCbCr.b(), 
//        //    yCbCr.r() + d * yCbCr.g()
//        //    );
//
//        return Color(
//                1.0f * yCbCr.r() + 0                            + 1.4746f * yCbCr.b(),
//                1.0f * yCbCr.r() - 0.16455312f * yCbCr.g()      - 0.57135312f * yCbCr.b(),
//                1.0f * yCbCr.r() + 0                            + 1.8814f * yCbCr.g()
//        );
//    }
//
//    /// @brief Given a Color rgb the pointers h,s,v contain the converted color.
//    void rgb2hsv(Color rgb, float* h, float* s, float* v)
//    {
//        float maxIntensity_ = max(max(rgb.r(), rgb.g()), rgb.b());
//        float minIntensity_ = min(min(rgb.r(), rgb.g()), rgb.b());
//
//        float delta = maxIntensity_ - minIntensity_;
//
//        unsigned int maxRGB = 0;
//        if(maxIntensity_ == rgb.r())
//            maxRGB = 1;
//        else if(maxIntensity_ == rgb.g())
//            maxRGB = 2;
//        else if(maxIntensity_ == rgb.b())
//            maxRGB = 3;
//        
//        switch (maxRGB)
//        {
//        case 1:
//            *h = 60 * fmod(((rgb.g() - rgb.b()) / delta), 6.0f);
//            break;
//        case 2:
//            *h = 60 * ((rgb.b() - rgb.r()) / delta + 2.0f);
//            break;
//        case 3:
//            *h = 60 * ((rgb.r() - rgb.g()) / delta + 4.0f);
//            break;
//        default:
//            break;
//        }
//
//        if(*h < 0)
//            *h = *h + 360.0f;
//
//        if(delta == 0)
//            *h = 0;
//
//        *s = maxIntensity_ == 0.0f ? 0.0f : delta / maxIntensity_;
//        *v = maxIntensity_;
//    }
//
//    /// @brief Given h,s,v the color pointer rgb contains the converted color.
//    void hsv2rgb(float h, float s, float v, Color* rgb)
//    {
//        float c = v * s;
//        float x = c * (1.0f - fabs(fmod(h / 60.0f, 2.0f) - 1.0f));
//        float m = v - c;
//
//        Color rgb_ = Color(-1);
//        if(h >= 0.0f && h < 60.0f)
//            rgb_ = Color(c, x, 0);
//        else if(h >= 60.0f && h < 120.0f)
//            rgb_ = Color(x, c, 0);
//        else if(h >= 120.0f && h < 180.0f)
//            rgb_ = Color(0, c, x);
//        else if(h >= 180.0f && h < 240.0f)
//            rgb_ = Color(0, x, c);
//        else if(h >= 240.0f && h < 300.0f)
//            rgb_ = Color(x, 0, c);
//        else if(h >= 300.0f && h < 360.0f)
//            rgb_ = Color(c, 0, x);
//
//        *rgb = Color(rgb_.r() + m, rgb_.g() + m, rgb_.b() + m);
//    }
//
//        /// @brief Given a Color rgb the pointers h,s,v contain the converted color.
//    void rgb2hsl(Color rgb, float* h, float* s, float* l)
//    {
//        float maxIntensity_ = max(max(rgb.r(), rgb.g()), rgb.b());
//        float minIntensity_ = min(min(rgb.r(), rgb.g()), rgb.b());
//        
//        float delta = maxIntensity_ - minIntensity_;
//
//        unsigned int maxRGB = 0;
//        if(maxIntensity_ == rgb.r())
//            maxRGB = 1;
//        else if(maxIntensity_ == rgb.g())
//            maxRGB = 2;
//        else if(maxIntensity_ == rgb.b())
//            maxRGB = 3;
//        
//        switch (maxRGB)
//        {
//        case 1:
//            *h = 60 * fmod(((rgb.g() - rgb.b()) / delta), 6.0f);
//            break;
//        case 2:
//            *h = 60 * ((rgb.b() - rgb.r()) / delta + 2.0f);
//            break;
//        case 3:
//            *h = 60 * ((rgb.r() - rgb.g()) / delta + 4.0f);
//            break;
//        default:
//            break;
//        }
//
//        if(*h < 0)
//            *h = *h + 360.0f;
//
//        if(delta == 0)
//            *h = 0;
//
//        
//
//        *l = (maxIntensity_ + minIntensity_) / 2;
//
//        //if(maxIntensity_ > 1.0f)
//        //{
//        //    *l /= maxIntensity_;
//        //    delta /= maxIntensity_;
//        //}
//
//        *s = delta == 0.0f ? 0.0f : *l == 1.0f ? delta / abs(1.0f - abs(2.0f * (*l - 0.001f) - 1.0f)) : delta / abs(1.0f - abs(2.0f * *l - 1.0f));
//
//        //if(maxIntensity_ > 1.0f)
//        //    *l *= maxIntensity_;
//    }
//
//    /// @brief Given h,s,v the color pointer rgb contains the converted color.
//    void hsl2rgb(float h, float s, float l, Color* rgb)
//    {
//        float c = l == 1.0f ? abs(1.0f - abs(2.0f * (l - 0.001f) - 1.0f)) * s : abs(1.0f - abs(2.0f * l - 1.0f)) * s;
//        float x = c * (1.0f - fabs(fmod(h / 60.0f, 2.0f) - 1.0f));
//        float m = abs(l - c / 2.0f);
//
//        Color rgb_ = Color(-1);
//        if(h >= 0.0f && h < 60.0f)
//            rgb_ = Color(c, x, 0);
//        else if(h >= 60.0f && h < 120.0f)
//            rgb_ = Color(x, c, 0);
//        else if(h >= 120.0f && h < 180.0f)
//            rgb_ = Color(0, c, x);
//        else if(h >= 180.0f && h < 240.0f)
//            rgb_ = Color(0, x, c);
//        else if(h >= 240.0f && h < 300.0f)
//            rgb_ = Color(x, 0, c);
//        else if(h >= 300.0f && h < 360.0f)
//            rgb_ = Color(c, 0, x);
//
//        *rgb = Color(rgb_.r() + m, rgb_.g() + m, rgb_.b() + m);
//        //*rgb = Color(m);
//    }