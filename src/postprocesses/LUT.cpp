#include <cstring>
#include <fstream>
#include <lightwave.hpp>

namespace lightwave {

struct LUT {
    int size;
    std::vector<Color> colors;
    //TODO set path for LUT

    float alex_strtof(char *head, char **endPtr) {
        while (*head == ' ')
            ++head;

        // bool isNegative = *head == '-';
        // if (isNegative) head++;

        /// @todo perhaps switch to int?
        long decimal      = 0;
        int baseExp       = 0;
        bool hasSeenPoint = false;

        while (true) {
            char chr = *(head++);
            if (chr >= '0' && chr <= '9') {
                /// @todo handle overflows of 'decimal' gracefully, maybe by
                /// limiting number of digits read
                decimal *= 10;
                decimal += chr - '0';
                baseExp -= hasSeenPoint ? 1 : 0;
            } else if (chr == '.') {
                assert(!hasSeenPoint);
                hasSeenPoint = true;
            } else if (chr == ' ' || chr == '\n') {
                break;
            } else {
                assert(!"invalid character");
            }
        }

        *endPtr = head - 1;

        /// @todo is this precise enough?
        // if (isNegative) decimal *= -1;
        auto baseExp2 = baseExp * baseExp;
        auto baseExp4 = baseExp2 * baseExp2;
        auto baseExp8 = baseExp4 * baseExp4;
        return 1.f * decimal * baseExp2 * baseExp8;
    }

    void readLUT(std::string filename) {
        std::ifstream in(filename);
        std::string contents((std::istreambuf_iterator<char>(in)),
                             std::istreambuf_iterator<char>());

        char *file = new char[contents.size() + 1];
        std::strcpy(file, contents.c_str());

        char *head    = &file[0];
        char **endPtr = &head;

        while (true) {
            if (*head == 'L' && (*(head - 1) == '\n' || head == &file[0])) {
                head += 12;
                size = static_cast<int>(alex_strtof(head, endPtr));
                colors.resize(size * size * size);
                head = *endPtr;

                break;
            } else {
                head++;
            }
        }

        while (!(*head >= '0' && *head <= '9')) {
            head++;
        }

        int x        = 0;
        int y        = 0;
        int z        = 0;
        int size_lut = size * size * size;
        for (int i = 0; i < size_lut; i++) {
            while (!(*head >= '0' && *head <= '9')) {
                head++;
            }

            if (x >= size) {
                x = 0;
                y++;
            }
            if (y >= size) {
                y = 0;
                z++;
            }

            Color color = Color(0);

            color.r() = alex_strtof(head, endPtr);
            color.g() = alex_strtof(head, endPtr);
            color.b() = alex_strtof(head, endPtr);

            colors[x * size * size + y * size + z] = color;
            x++;
        }

        delete file;
    }

    const Color &get(const int x, const int y, const int z) const {
        return colors[x * size * size + y * size + z];
    }

    Color trilinear(Color rgb) {
        rgb *= static_cast<float>(size - 1);

        int rl = int(floor(rgb.r()));
        int gl = int(floor(rgb.g()));
        int bl = int(floor(rgb.b()));
        int ru = int(ceil(rgb.r()));
        int gu = int(ceil(rgb.g()));
        int bu = int(ceil(rgb.b()));

        Color c0 = get(rl, gl, bl);
        Color c1 = get(rl, gl, bu);
        Color c2 = get(rl, gu, bl);
        Color c3 = get(rl, gu, bu);
        Color c4 = get(ru, gl, bl);
        Color c5 = get(ru, gl, bu);
        Color c6 = get(ru, gu, bl);
        Color c7 = get(ru, gu, bu);

        Color l1   = lerp(c0, c1, rgb.b() - bl);
        Color l2   = lerp(c2, c3, rgb.b() - bl);
        Color l3   = lerp(c4, c5, rgb.b() - bl);
        Color l4   = lerp(c6, c7, rgb.b() - bl);
        Color ll12 = lerp(l1, l2, rgb.g() - gl);
        Color ll34 = lerp(l3, l4, rgb.g() - gl);
        return lerp(ll12, ll34, rgb.r() - rl);
    }
};

class ColorLUTs : public Postprocess {
    LUT lut;

public:
    ColorLUTs(const Properties &properties) : Postprocess(properties) {
        lut.readLUT("../LUTs/Base sRGB.cube");
    }

    /// @brief Image gets more/less contrast/brightness
    void execute() override {
        // Output initialized to resolution of input
        m_output->initialize(m_input->resolution());

        // Loop over each Pixel
        for (int y = 0; y < m_output->resolution().y(); y++) {
            for (int x = 0; x < m_output->resolution().x(); x++) {
                Point2i pixel(x, y);
                Color rgb = m_input->get(pixel);

                for (int i = 0; i < 3; i++) {
                    // rgb[i] = log10f(rgb[i]+1);
                    // rgb[i] = powf(rgb[i], 2.4f);
                }

                // check wether Color is above 1.0f
                float maxIntensity = 0.0f;
                bool isEmissive    = false;
                if (rgb.r() > 1.0f || rgb.g() > 1.0f || rgb.b() > 1.0f) {
                    maxIntensity = max(max(rgb.r(), rgb.g()), rgb.b());
                    isEmissive   = true;
                }

                if (isEmissive)
                    rgb /= (maxIntensity + 0.001f);

                rgb = lut.trilinear(rgb);

                if (isEmissive)
                    rgb *= (maxIntensity + 0.001f);

                for (int i = 0; i < 3; i++) {
                    // rgb[i] = log10f(rgb[i]+1);
                    // rgb[i] = powf(rgb[i], 1.0f/2.4f);
                }

                m_output->get(pixel) = rgb;
            }
        }

        m_output->save();
    }

    std::string toString() const override {
        return tfm::format(
            "LUT[\n"
            "  input = %s,\n"
            "  output = %s,\n"
            "]",
            indent(m_input),
            indent(m_output)
        );
    }
};

} // namespace lightwave

REGISTER_POSTPROCESS(ColorLUTs, "LUT");
