#include <lightwave.hpp>

namespace lightwave {

enum InterpolType { EASE, CARDINAL, LINEAR, B_SPLINE, CONSTANT };
enum ColorSpace { RGB, HSV, HSL };
enum HueInterpol { NEAR, FAR, CW, CCW };

class ColorRamp : public Texture {

public:
    class Element : public Object {

    protected:
        Color m_color;
        float m_alpha;
        float m_position;

    public:
        Element() {
            m_color    = Color(0);
            m_position = 0;
            m_alpha    = 1;
        }

        Element(const Color color, const float position, float alpha = 1) {
            m_color    = color;
            m_position = position;
            m_alpha    = alpha;
        }

        Element(const Properties &properties) {
            m_color    = properties.get<Color>("color");
            m_position = properties.get<float>("position");
            m_alpha    = properties.get<float>("alpha");
        }

        float position() const { return m_position; }
        Color color() const { return m_color; }
        float alpha() const { return m_alpha; }

        std::string toString() const override {
            return tfm::format(
                "Element[\n"
                "color = %s\n",
                "position = %s\n",
                "alpha = %s\n",
                "]",
                indent(m_color),
                indent(m_position),
                indent(m_alpha));
        }
    };

protected:
    ref<Texture> m_factor;
    ColorSpace m_colorMode;
    InterpolType m_interpolationType;
    HueInterpol m_hueInterpolation;
    std::vector<Element> m_elements = {};
    int numElems;

public:
    ColorRamp(const Properties &properties) {

        m_factor = properties.get<Texture>("factor");

        m_colorMode = properties.getEnum<ColorSpace>(
            "colorMode", { { "RGB", RGB }, { "HSV", HSV }, { "HSL", HSL } });
        m_interpolationType =
            properties.getEnum<InterpolType>("interpolationType",
                                             { { "EASE", EASE },
                                               { "CARDINAL", CARDINAL },
                                               { "LINEAR", LINEAR },
                                               { "B_SPLINE", B_SPLINE },
                                               { "CONSTANT", CONSTANT } });
        m_hueInterpolation = properties.getEnum<HueInterpol>(
            "hueInterpolation",
            { { "NEAR", NEAR }, { "FAR", FAR }, { "CW", CW }, { "CCW", CCW } });

        auto elements = properties.getChildren<Element>();
        numElems      = elements.size();
        assert_condition(numElems > 0,
                         logger(EError, "ColorRamp needs to have elements!"););

        for (int i = 0; i < numElems; i++) {
            m_elements.push_back(*elements[i]);
        }
        std::sort(m_elements.begin(),
                  m_elements.end(),
                  [](const Element &el1, const Element &el2) {
                      return el1.position() < el2.position();
                  });
    }

    Color evaluate(const Point2 &uv, const Context &cont) const override {

        assert_condition(
            m_colorMode == RGB,
            logger(EError, "ColorRamp only support color mode RGB!"););

        float fac = clamp(m_factor->scalar(uv, cont), 0.f, 1.f);

        if (m_interpolationType == CONSTANT || m_interpolationType == LINEAR ||
            m_interpolationType == EASE) {
            if (m_elements[0].position() >= fac)
                return m_elements[0].color();
            if (m_elements[numElems - 1].position() <= fac)
                return m_elements[numElems - 1].color();
        }

        Element l, r;
        float t = 0;
        int i;
        if (m_elements[0].position() >= fac) { // We only enter the if parts if
                                               // we do CARDINAL or B_SPLINE
            i = 0;
            r = m_elements[0];
            l = Element(r.color(), 0.f);
        } else if (m_elements[numElems - 1].position() <= fac) {
            i = numElems;
            l = m_elements[numElems - 1];
            r = Element(l.color(), 1.f);
        } else {
            i = 1;
            while (i < numElems) {
                if (m_elements[i].position() >= fac) {
                    l = m_elements[i - 1];
                    r = m_elements[i];
                    break;
                }
                i++;
            }
        }
        if (l.position() != r.position())
            t = (fac - l.position()) / (r.position() - l.position());
        else
            t = (i == numElems) ? 1.f
                                : 0.f; // Only needed for CARDINAL and B_SPLINE

        switch (m_interpolationType) {
        case CONSTANT:
            return l.color();
        case EASE:
            return interpolateEase(t, l.color(), r.color());
        case LINEAR:
            return interpolateLinear(t, l.color(), r.color());
        default: { // We are in the cases CARDINAL or B_SPLINE and interpolate
                   // the 4 neighboring colors, see colorband.cc
            Element ll, rr; // LeftLeft and RightRight neighbors
            if (i >= numElems - 1)
                rr = r;
            else
                rr = m_elements[i + 1];

            if (i < 2)
                ll = l;
            else
                ll = m_elements[i - 1];

            const Colors4 colors = {
                rr.color(), r.color(), l.color(), ll.color()
            };
            if (m_interpolationType == B_SPLINE)
                return interpolateCubicBSpline(
                    t, colors); // Blender:
                                // key_curve_position_weights(fac,
                                // t, KEY_BSPLINE);

            const Floats4 points = toFloats4(
                { rr.position(), r.position(), l.position(), ll.position() });
            return interpolateCatmullRom(
                t, colors, points); // Blender: key_curve_position_weights(fac,
                                    // t, KEY_CARDINAL);
        }
        }
        // This should be unreachable!
        return Color(0);
    }

    std::string toString() const override {
        return tfm::format(
            "ColorRamp[\n"
            "  factor = %s\n"
            "  colorMode = %s\n"
            "  interpolationType = %s\n"
            "  hueInterpolation = %s\n"
            "]",
            indent(m_factor),
            indent(m_colorMode),
            indent(m_interpolationType),
            indent(m_hueInterpolation));
    }
};

} // namespace lightwave

REGISTER_TEXTURE(ColorRamp, "color_ramp")
REGISTER_CLASS(ColorRamp::Element, "color_ramp_element", "default")
