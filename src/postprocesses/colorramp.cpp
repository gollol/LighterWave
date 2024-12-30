#include <lightwave.hpp>

namespace lightwave {

enum InterpolType{EASE, CARDINAL, LINEAR, B_SPLINE, CONSTANT};
enum ColorSpace{RGB, HSV, HSL};
enum HueInterpol{NEAR, FAR, CW, CCW};

class Element : public Object {

protected:
    Color m_color;
    float m_alpha;
    float m_position;

public:
    Element(){
        m_color = Color(0);
        m_position = 0;
        m_alpha = 1;
    }

    Element(const Color color, const float position, float alpha = 1){
        m_color = color;
        m_position = position;
        m_alpha = alpha;
    }

    Element(const Properties &properties) {
        m_color = properties.get<Color>("color");
        m_position = properties.get<float>("position");
        m_alpha = properties.get<float>("alpha");
    }

    float position() const {return m_position;}
    Color color() const {return m_color;}
    float alpha() const {return m_alpha; }

    std::string toString() const override {
        return tfm::format("Element[\n"
                              "color = %s\n",
                              "position = %s\n",
                              "alpha = %s\n",
                           "]",
                           indent(m_color), indent(m_position), indent(m_alpha));
    }
};

class PostColorRamp : public Postprocess {

public:

protected:
    ref<Texture> m_factor;
    ColorSpace m_colorMode;
    InterpolType m_interpolationType;
    HueInterpol m_hueInterpolation;
    std::vector<Element> m_elements = {};
    int numElems;

public:
    PostColorRamp(const Properties &properties) : Postprocess(properties) 
    {
        m_factor = properties.get<Texture>("factor");

        m_colorMode = properties.getEnum<ColorSpace>("colorMode", {
            {"RGB", RGB}, {"HSV", HSV}, {"HSL", HSL}});
        m_interpolationType = properties.getEnum<InterpolType>("interpolationType", {
            {"EASE", EASE}, {"CARDINAL", CARDINAL}, {"LINEAR", LINEAR}, {"B_SPLINE", B_SPLINE}, {"CONSTANT", CONSTANT}});
        m_hueInterpolation = properties.getEnum<HueInterpol>("hueInterpolation", {
            {"NEAR", NEAR}, {"FAR", FAR}, {"CW", CW}, {"CCW", CCW}});

        auto elements = properties.getChildren<Element>();
        numElems = elements.size();
        assert(numElems>0 && "PostColorRamp needs to have elements!");

        for(int i = 0; i<numElems; i++){
            m_elements.push_back(*elements[i]);
        }
        std::sort(m_elements.begin(), m_elements.end(), [](const Element &el1, const Element &el2) {return el1.position() < el2.position();});
    }

    void execute() override 
    {
    //Color evaluate(const Point2 &uv, const Context &cont) const override {
        m_output->initialize(m_input->resolution());


        if(m_colorMode != RGB){
            lightwave_throw("Only RGB color format supported in PostColorRamp so far!");
            return;
        }

        for(int y = 0; y < m_output->resolution().y(); y++)
        for(int x = 0; x < m_output->resolution().x(); x++)
        {
            Point2i pixel(x,y);
            Point2 uv = Point2((float)x/m_output->resolution().x(), (float)(m_output->resolution().y() - y) / m_output->resolution().y());
            Color out = Color(0);

            //float fac = clamp(m_factor->scalar(uv), 0.f, 1.f);
            float fac = m_input->get(pixel).r();

            if(m_interpolationType == CONSTANT || m_interpolationType == LINEAR || m_interpolationType == EASE)
            {
                if(m_elements[0].position() >= fac)
                    out = m_elements[0].color();
                if(m_elements[numElems-1].position() <= fac)
                    out = m_elements[numElems-1].color();
            }

            Element l, r;
            float t = 0;
            int i;
            if(m_elements[0].position() >= fac){        //We only enter the if parts if we do CARDINAL or B_SPLINE
                i = 0;
                r = m_elements[0];
                l = Element(r.color(), 0.f);
            }else if(m_elements[numElems-1].position() <= fac){
                i = numElems;
                l = m_elements[numElems-1];
                r = Element(l.color(), 1.f);
            }else{
                i = 1;
                while(i<numElems){
                    if(m_elements[i].position() >= fac){
                        l = m_elements[i-1];
                        r = m_elements[i];
                        break;
                    }
                    i++;
                }
            }

            if(l.position()!=r.position())
                t = (fac - l.position())/(r.position() - l.position());
            else
                t = (i==numElems) ? 1.f : 0.f;    // Only needed for CARDINAL and B_SPLINE

            switch(m_interpolationType){
                case CONSTANT:
                    out = l.color();
                    break;
                case EASE:
                    t = t * t * (3.f - 2.f * t);    //We don't return because we want to interpolate linearly now!
                case LINEAR:
                    out = (1-t) * l.color() + t * r.color();
                    break;
                default: {  //We are in the cases CARDINAL or B_SPLINE and interpolate the 4 neighboring colors, see colorband.cc
                    //Element ll, rr;     //LeftLeft and RightRight neighbors
                    //if (i>=numElems-1) rr = r;
                    //else               rr = m_elements[i+1];
//
                    //if (i<2) ll = l;
                    //else     ll = m_elements[i-1];
//
                    //const Colors4 colors = {rr.color(), r.color(), l.color(), ll.color()};
                    //if(m_interpolationType == B_SPLINE)
                    //    return interpolateCubicBSpline(t, colors);    //Blender: key_curve_position_weights(fac, t, KEY_BSPLINE);
//
                    //const Floats4 points = toFloats4({rr.position(), r.position(), l.position(), ll.position()});
                    //return interpolateCatmullRom(t, colors, points);  //Blender: key_curve_position_weights(fac, t, KEY_CARDINAL);
                }
            }
            // This should be unreachable!

            m_output->get(pixel) = out;
        }

        m_output->save();
    }

    std::string toString() const override {
        return tfm::format("PostColorRamp[\n"
                           "  m_factor = %s\n"
                           "  m_colorMode = %s\n"
                           "  m_interpolationType = %s\n"
                           "  m_hueInterpolation = %s\n"
                           "]",
                           indent(m_factor),
                           indent(m_colorMode),
                           indent(m_interpolationType),
                           indent(m_hueInterpolation));
    }
};

} // namespace lightwave

REGISTER_POSTPROCESS(PostColorRamp, "colorramp")
REGISTER_CLASS(Element, "colorRampElement", "default")
