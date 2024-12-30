#include <lightwave.hpp>

namespace lightwave {

class MipImageTexture : public Texture {
    enum class BorderMode {
        Clamp,
        Repeat,
    };

    enum class FilterMode {
        Nearest,
        Bilinear,
    };

    std::vector<ref<Image>> m_images;
    float m_exposure;
    BorderMode m_border;
    FilterMode m_filter;
    int m_mip;

    void buildMipmaps() {
        int levels = static_cast<uint32_t>(std::floor(
                         std::log2(std::max(m_images[0]->resolution().x(),
                                            m_images[0]->resolution().y())))) +
                     1;

        int mipWidth  = m_images[0]->resolution().x();
        int mipHeight = m_images[0]->resolution().y();
        for (int i = 1; i < levels; i++) {
            m_images.push_back(std::make_shared<Image>(
                Point2i(mipWidth > 1 ? mipWidth / 2 : 1,
                        mipHeight > 1 ? mipHeight / 2 : 1)));

            for (int y = 0; y < (mipHeight > 1 ? mipHeight / 2 : 1); y++) {
                for (int x = 0; x < (mipWidth > 1 ? mipWidth / 2 : 1); x++) {
                    int x1 = (x * 2) + 0;
                    int y1 = (y * 2) + 0;
                    int x2 = x1 + 1 >= m_images[i - 1]->resolution().x()
                                 ? x1
                                 : x1 + 1;
                    int y2 = y1 + 1 >= m_images[i - 1]->resolution().y()
                                 ? y1
                                 : y1 + 1;

                    Color a = m_images[i - 1]->get(Point2i(x1, y1));
                    Color b = m_images[i - 1]->get(Point2i(x2, y1));
                    Color c = m_images[i - 1]->get(Point2i(x1, y2));
                    Color d = m_images[i - 1]->get(Point2i(x2, y2));

                    m_images[i]->get(Point2i(x, y)) = (a + b + c + d) / 4.0f;
                }
            }

            mipWidth /= 2;
            mipHeight /= 2;
        }

#if 0
            for (int i = 0; i < m_images.size(); i++) {
                m_images[i]->saveAt("height_mipmap_" + std::to_string(i) + ".exr");
            }
#endif
    }

    static inline int clampBorder(int x, int w) {
        return std::clamp(x, 0, w - 1);
    }

    static inline int repeatBorder(int x, int w) {
        const int t = x % w;
        return t < 0 ? t + w : t;
    }

    inline Point2i handleBorder(const Point2i &uv, int mip) const {
        const auto res = m_images[mip]->resolution();
        switch (m_border) {
        default:
        case BorderMode::Clamp:
            return Point2i(clampBorder(uv.x(), res.x()),
                           clampBorder(uv.y(), res.y()));
        case BorderMode::Repeat:
            return Point2i(repeatBorder(uv.x(), res.x()),
                           repeatBorder(uv.y(), res.y()));
        }
    }

    struct ImageTexelPoint {
        Point2i i;
        Point2 f;
    };

    inline ImageTexelPoint mapUV(const Point2 &uv, int mip,
                                 const Vector2 &offset = Vector2(0.5f)) const {
        const auto res = Vector2i(m_images[mip]->resolution());
        const Vector2 k =
            Vector2(uv.x(), 1 - uv.y()) * res.cast<float>() - offset;
        const Point2i i =
            Point2i((int) std::floor(k.x()), (int) std::floor(k.y()));
        const Point2 f = k - i.cast<float>();
        return ImageTexelPoint{ i, f };
    }

    inline Color nearestFilter(const Point2 &uv, int mip) const {
        const auto texel =
            mapUV(uv, mip, Vector2(0)); // Nearest filter uses no offset
        // (Cycles/OpenImageIO convention)
        const auto p = handleBorder(texel.i, mip);
        return (*m_images[mip])(p);
    }

    inline Color bilinearFilter(const Point2 &uv, int mip) const {
        const auto texel = mapUV(uv, mip);
        const auto p0    = handleBorder(texel.i, mip);
        const auto p1    = handleBorder(Vector2i(texel.i) + Vector2i(1), mip);

        const auto v00 = (*m_images[mip])(p0);
        const auto v01 = (*m_images[mip])(Point2i(p0.x(), p1.y()));
        const auto v10 = (*m_images[mip])(Point2i(p1.x(), p0.y()));
        const auto v11 = (*m_images[mip])(p1);

        return lerp(lerp(v00, v10, texel.f.x()),
                    lerp(v01, v11, texel.f.x()),
                    texel.f.y());
    }

public:
    MipImageTexture(const Properties &properties) {
        if (properties.has("filename")) {
            m_images.push_back(std::make_shared<Image>(properties));
        } else {
            m_images.push_back(properties.getChild<Image>());
        }
        buildMipmaps();
        m_exposure = properties.get<float>("exposure", 1);
        m_mip      = properties.get<int>("mip", 0);
        assert_condition(m_mip >= 0, {});
        assert_condition(m_mip < (int) m_images.size(), {});

        // clang-format off
            m_border = properties.getEnum<BorderMode>("border", BorderMode::Repeat, {
                { "clamp", BorderMode::Clamp },
                { "repeat", BorderMode::Repeat },
                });

            m_filter = properties.getEnum<FilterMode>("filter", FilterMode::Bilinear, {
                { "nearest", FilterMode::Nearest },
                { "bilinear", FilterMode::Bilinear },
                });
        // clang-format on
    }

    Color evaluate(const Point2 &uv, const Context &cont) const override {
        return evaluate(uv, cont, m_mip);
    }

    Color evaluate(const Point2 &uv, const Context &cont,
                   int mip) const override {
        assert_condition(mip >= 0, {});
        assert_condition(mip < (int) m_images.size(), {});

        switch (m_filter) {
        default:
        case FilterMode::Bilinear:
            return m_exposure * bilinearFilter(uv, mip);
        case FilterMode::Nearest:
            return m_exposure * nearestFilter(uv, mip);
        }
    }

    std::string toString() const override {
        return tfm::format(
            "ImageTexture[\n"
            "  image = %s,\n"
            "  mips = %s,\n"
            "  exposure = %f,\n"
            "]",
            indent(m_images[0]),
            indent(m_images.size()),
            m_exposure);
    }
};

} // namespace lightwave

REGISTER_TEXTURE(MipImageTexture, "mipImage")
