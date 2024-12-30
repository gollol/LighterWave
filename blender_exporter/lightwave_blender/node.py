import bpy
import os

from .utils import *
from .addon_preferences import get_prefs
from .xml_node import XMLNode
from .registry import SceneRegistry, Token
from .export_context import ExportContext
from .node_graph import RMInput, RMNode


def _export_fallback():
    return XMLNode("texture", type="constant", value=0)

def _export_default(export_ctx: ExportContext, input: RMInput, exposure):
    if not input.has_value():
        return _export_fallback()

    factor = exposure or 1
    value = input.value
    if isinstance(value, float):
        value *= factor
    elif isinstance(value, list):
        value = [ v * factor for v in value[0:3] ]
    return XMLNode("texture", type="constant", value=str_flat_array(value))

def _export_scalar_value(export_ctx: ExportContext, input: RMInput, exposure):
    factor = exposure or 1
    return XMLNode("texture", type="constant", value=input.linked_node().bl_node.outputs[0].default_value * factor)

def _export_rgb_value(export_ctx: ExportContext, input: RMInput, exposure):
    factor = exposure or 1
    node = input.linked_node().bl_node
    default_value = [ v * factor for v in node.outputs[0].default_value[0:3] ]
    return XMLNode("texture", type="constant", value=str_flat_array(default_value))

def _export_image(export_ctx: ExportContext, image, path, is_f32=False, keep_format=False):
    if os.path.exists(path) and not export_ctx.settings.overwrite_existing_textures:
        return

    # This should not be needed, maybe??
    ## Make sure the image is loaded to memory, so we can write it out
    #if not image.has_data:
    #    image.pixels[0]

    # Export the actual image data
    try:
        old_path = image.filepath_raw
        old_format = image.file_format
        try:
            image.filepath_raw = path
            if not keep_format:
                image.file_format = "PNG" if not is_f32 else "OPEN_EXR"
            image.save()
        finally:  # Never break the scene!
            image.filepath_raw = old_path
            image.file_format = old_format
    except:
        if not keep_format:
            raise RuntimeError(
                "Can not change the original format of given image")
        # Try other way
        image.save_render(path)


def _handle_image(export_ctx: ExportContext, image: bpy.types.Image):
    tex_dir_name = get_prefs().tex_dir_name

    if image.source == 'GENERATED':
        img_name = image.name + \
            (".png" if not image.use_generated_float else ".exr")
        img_path = os.path.join(tex_dir_name, img_name)
        _export_image(export_ctx, image, os.path.join(export_ctx.path, img_path), is_f32=image.use_generated_float)
        return img_path.replace('\\', '/') # Ensure the image path is not using \ to keep the xml valid
    elif image.source == 'FILE':
        filepath = image.filepath_raw if image.filepath_raw is not None else image.filepath
        img_path = bpy.path.abspath(bpy.path.resolve_ncase(filepath), library=image.library)
        try:
            img_path = bpy.path.relpath(img_path, start=export_ctx.path)
        except:
            # relpath can fail on Windows if paths have different drive letters
            print("unable to create relative path")
        img_path = img_path.replace("\\", "/")
        if img_path.startswith("//"):
            img_path = img_path[2:]

        copy_image = image.packed_file or getattr(
            export_ctx.settings, "copy_images", False) or img_path == ''

        if copy_image:
            export_ctx.debug(f"Copying image {img_path}") 
            img_name = bpy.path.basename(img_path)

            # Special case: We can not export PNG if bit depth is not 8 (or 32), for whatever reason
            if img_name == '' or image.depth > 32 or image.depth == 16:
                keep_format = False
                if image.depth > 32 or image.depth == 16 or image.file_format in ["OPEN_EXR", "OPEN_EXR_MULTILAYER", "HDR"]:
                    is_f32 = True
                    extension = ".exr"
                else:
                    is_f32 = False
                    extension = ".png"
                img_path = os.path.join(tex_dir_name, image.name + extension)
            else:
                keep_format = True
                is_f32 = False  # Does not matter
                img_path = os.path.join(tex_dir_name, img_name)

            try:
                _export_image(export_ctx, image, os.path.join(export_ctx.path, img_path),
                                is_f32=is_f32, keep_format=keep_format)
            except:
                # Above failed, so give this a try
                img_path = os.path.join(tex_dir_name, img_name)
                _export_image(export_ctx, image, os.path.join(export_ctx.path, img_path),
                                  is_f32=False, keep_format=True)
        return img_path.replace('\\', '/') # Ensure the image path is not using \ to keep the xml valid
    else:
        export_ctx.error(f"Image type {image.source} not supported")
        return None


def _export_image_texture(export_ctx: ExportContext, input: RMInput, exposure):
    bl_node: bpy.types.ShaderNodeTexImage = input.linked_node().bl_node
    if not bl_node.image:
        export_ctx.error(f"Image node {bl_node.name} has no image")
        return _export_fallback()

    def export():
        img_path = _handle_image(export_ctx, bl_node.image)
        kw = {}

        if exposure is not None and exposure != 1:
            kw["exposure"] = exposure
        
        if bl_node.extension == "EXTEND" or bl_node.extension == "CLIP":
            kw["border"] = "clamp"

        if bl_node.interpolation == "Closest":
            kw["filter"] = "nearest"

        cs = bl_node.image.colorspace_settings.name
        if cs == "Linear" or cs == "Non-Color" or cs == "Raw":
            kw["linear"] = True
        elif cs == "sRGB":
            pass
        else:
            export_ctx.error(f"Unsupported color space {cs}")

        return XMLNode("texture", type="mipImage", filename=img_path, **kw)
    
    return export_ctx.registry.export(bl_node.image.original, export)


def _export_env_image_texture(export_ctx: ExportContext, input: RMInput, exposure):
    # Handle it the same way as standard textures, but with tiny differences (e.g., fixed border handling)
    bl_node: bpy.types.ShaderNodeTexEnvironment = input.linked_node().bl_node
    if not bl_node.image:
        export_ctx.error(f"Image node {bl_node.name} has no image")
        return _export_fallback()

    def export():
        img_path = _handle_image(export_ctx, bl_node.image)
        kw = {}

        if exposure is not None and exposure != 1:
            kw["exposure"] = exposure
        
        kw["border"] = "clamp"

        if bl_node.interpolation == "Closest":
            kw["filter"] = "nearest"

        cs = bl_node.image.colorspace_settings.name
        if cs == "Linear" or cs == "Non-Color" or cs == "Raw":
            kw["linear"] = True
        elif cs == "sRGB":
            pass
        else:
            export_ctx.error(f"Unsupported color space {cs}")

        return XMLNode("texture", type="image", filename=img_path, **kw)
    
    return export_ctx.registry.export(bl_node.image.original, export)

def _export_normal_map(export_ctx: ExportContext, input: RMInput, exposure, **args):
    return export_node(export_ctx, input.linked_node().input("Color"), exposure, **args) # TODO: Strength


def _export_val_to_rgb(export_ctx: ExportContext, input: RMInput, exposure):
    factor = exposure or 1
    node: RMNode = input.linked_node()
    output: str = input.linked_output()
    bl_node: bpy.types.ShaderNodeValToRGB = node.bl_node

    ramp = bl_node.color_ramp

    xml_node = XMLNode("texture", type="color_ramp", colorMode=ramp.color_mode, interpolationType=ramp.interpolation, hueInterpolation=ramp.hue_interpolation)
    xml_node.add_child(export_node(export_ctx, node.input("Fac")), name="factor")

    for element in ramp.elements.values():
        element_xml_node = XMLNode("color_ramp_element", alpha=element.alpha, position=element.position)
        element_xml_node.add_child(XMLNode("color", value=str_flat_array([v * (exposure or 1) for v in element.color][:3])), name="color")

        xml_node.add_child(element_xml_node)
    
    return xml_node

def _export_curves_info(export_ctx: ExportContext, input: RMInput, unknown):
    var = input.linked_output().lower()
    if input.linked_output() == "Is Strand":
        var = "isStrand"
    elif input.linked_output() == "Tangent Normal":
        var = "tangentNormal"
    node = XMLNode("texture", type="curves_info", variable=var)
    return node

def _export_object_info(export_ctx: ExportContext, input: RMInput, unknown):
    var = input.linked_output().lower()
    if input.linked_output() == "Object Index":
        var = "objectIndex"
    elif input.linked_output() == "Material Index":
        var = "materialIndex"
    node = XMLNode("texture", type ="object_info", variable = var)
    return node

def _export_particle_info(export_ctx: ExportContext, input: RMInput, unknown):
    var = input.linked_output().lower()
    if input.linked_output() == "Angular Velocity":
        var = "angularVelocity"
    node = XMLNode("texture", type ="particle_info", variable = var)
    return node

def _export_geometry_info(export_ctx: ExportContext, input: RMInput, unknown):
    var = input.linked_output().lower()
    if input.linked_output() == "True Normal":
        var = "trueNormal"
    elif input.linked_output() == "Random Per Island":
        var = "randomPerIsland"
    node = XMLNode("texture", type ="geometry_info", variable = var)
    return node

def _export_checkerboard_texture(export_ctx: ExportContext, input: RMInput, exposure):
    factor = exposure or 1
    node = input.linked_node()

    xml_node = XMLNode("texture", type="checkerboard")

    xml_node.add_child(export_node(export_ctx, node.input("Color1")), name="color0")
    xml_node.add_child(export_node(export_ctx, node.input("Color2")), name="color1")
    xml_node.add_child(export_node(export_ctx, node.input("Scale")), name="scale")

    return xml_node

def _export_hue_saturation(export_ctx: ExportContext, input: RMInput, exposure):
    factor = exposure or 1
    node = input.linked_node()

    xml_node = XMLNode("texture", type="hue_saturation")

    xml_node.add_child(export_node(export_ctx, node.input("Hue")), name="hue")
    xml_node.add_child(export_node(export_ctx, node.input("Saturation")), name="saturation")
    xml_node.add_child(export_node(export_ctx, node.input("Value")), name="value")
    xml_node.add_child(export_node(export_ctx, node.input("Fac")), name="factor")
    xml_node.add_child(export_node(export_ctx, node.input("Color")))

    return xml_node

def _export_invert(export_ctx: ExportContext, input: RMInput, exposure):
    node = input.linked_node()

    xml_node = XMLNode("texture", type="invert")
    xml_node.add_child(export_node(export_ctx, node.input("Fac")), name="fac")
    xml_node.add_child(export_node(export_ctx, node.input("Color")))

    return xml_node

def _export_combine_xyz(export_ctx: ExportContext, input: RMInput, exposure):
    node = input.linked_node()

    xml_node = XMLNode("texture", type="combine_xyz")
    xml_node.add_child(export_node(export_ctx, node.input("X")), name="x")
    xml_node.add_child(export_node(export_ctx, node.input("Y")), name="y")
    xml_node.add_child(export_node(export_ctx, node.input("Z")), name="z")

    return xml_node

def _export_fresnel(export_ctx: ExportContext, input: RMInput, exposure):
    node = input.linked_node()

    xml_node = XMLNode("texture", type="fresnel")
    xml_node.add_child(export_node(export_ctx, node.input("IOR")), name="ior")

    return xml_node

def _export_combine_xyz(export_ctx: SceneRegistry, input: RMInput, exposure):
    node = input.linked_node()

    xml_node = XMLNode("texture", type="combine_xyz")
    xml_node.add_child(export_node(export_ctx, node.input("X")), name="x")
    xml_node.add_child(export_node(export_ctx, node.input("Y")), name="y")
    xml_node.add_child(export_node(export_ctx, node.input("Z")), name="z")

    return xml_node

def _export_fresnel(export_ctx: ExportContext, input: RMInput, exposure):
    node = input.linked_node()

    xml_node = XMLNode("texture", type="fresnel")
    xml_node.add_child(export_node(export_ctx, node.input("IOR")), name="ior")

    return xml_node

def _export_brightness_contrast(export_ctx: ExportContext, input: RMInput, exposure):
    node = input.linked_node()

    xml_node = XMLNode("texture", type="brightnesscontrast")
    xml_node.add_child(export_node(export_ctx, node.input("Bright")), name="bright")
    xml_node.add_child(export_node(export_ctx, node.input("Contrast")), name="contrast")
    xml_node.add_child(export_node(export_ctx, node.input("Color")))

    return xml_node

supported_math_operations = {"ADD", "SUBTRACT", "MULTIPLY", "MULTIPLY_ADD", "MINIMUM", "MAXIMUM", "POWER", "GREATER_THAN", "LESS_THAN"}
def _export_math(export_ctx: ExportContext, input: RMInput, exposure):
    node = input.linked_node()
    blnode = node.bl_node
    operation = blnode.operation

    if operation not in supported_math_operations:
        export_ctx.warn("For the math node the " + str(operation) + " operation is not supported!")

    xml_node = XMLNode("texture", type="math", operation=operation, clamp=blnode.use_clamp)
    xml_node.add_child(export_node(export_ctx, node.input("Value")), name="value")
    xml_node.add_child(export_node(export_ctx, node.input("Value_001")), name="value_001")
    xml_node.add_child(export_node(export_ctx, node.input("Value_002")), name="value_002")

    return xml_node

supported_map_operations = {"LINEAR"}
supported_map_types = {"FLOAT"}
def _export_map_range(export_ctx: ExportContext, input: RMInput, exposure):
    node = input.linked_node()
    blnode = node.bl_node
    interpolation_type = blnode.interpolation_type
    data_type = blnode.data_type

    if interpolation_type not in supported_map_operations:
        export_ctx.warn("For the map_range node the " + str(interpolation_type) + " operation is not supported!")
    elif data_type not in supported_map_types:
        export_ctx.warn("For the map_range node the " + str(data_type) + " data type is not supported!")

    xml_node = XMLNode("texture", type="map_range", data_type=data_type, interpolation_type=interpolation_type, clamp=blnode.clamp)
    xml_node.add_child(export_node(export_ctx, node.input("Value")), name="value")
    xml_node.add_child(export_node(export_ctx, node.input("From Min")), name="fromMin")
    xml_node.add_child(export_node(export_ctx, node.input("From Max")), name="fromMax")
    xml_node.add_child(export_node(export_ctx, node.input("To Min")), name="toMin")
    xml_node.add_child(export_node(export_ctx, node.input("To Max")), name="toMax")

    return xml_node

supported_blend_modes = {"MIX", "MULTIPLY", "ADD"}
def _export_mix(export_ctx: ExportContext, input: RMInput, exposure):
    node = input.linked_node()
    blnode = node.bl_node

    if blnode.data_type == "RGBA" and blnode.blend_type not in supported_blend_modes:
        export_ctx.warn("For the mix node the " + str(blnode.blend_type) + " blend mode is not supported!")

    xml_node = XMLNode("texture", type="mix", data_type=blnode.data_type, factor_mode=blnode.factor_mode, clamp_factor=blnode.clamp_factor, clamp_result=blnode.clamp_result, blend_type=blnode.blend_type)
    xml_node.add_child(export_node(export_ctx, node.input("Factor_Float")), name="factor_Float")
    xml_node.add_child(export_node(export_ctx, node.input("Factor_Vector")), name="factor_Vector")
    xml_node.add_child(export_node(export_ctx, node.input("A_Color")), name="a_Color")
    xml_node.add_child(export_node(export_ctx, node.input("B_Color")), name="b_Color")
    xml_node.add_child(export_node(export_ctx, node.input("A_Float")), name="a_Float")
    xml_node.add_child(export_node(export_ctx, node.input("B_Float")), name="b_Float")
    xml_node.add_child(export_node(export_ctx, node.input("A_Vector")), name="a_Vector")
    xml_node.add_child(export_node(export_ctx, node.input("B_Vector")), name="b_Vector")

    return xml_node

def _export_noise_texture(export_ctx: ExportContext, input: RMInput, exposure):
    factor = exposure or 1
    node = input.linked_node()
    output = input.linked_output()
    blnode = node.bl_node

    uses_color = (output == "Color")

    xml_node = XMLNode("texture", type="noise", noisetype=blnode.noise_type, normalize=blnode.normalize, colorful=uses_color, dimensions=blnode.noise_dimensions)

    xml_node.add_child(export_node(export_ctx, node.input("Scale")), name="scale")
    xml_node.add_child(export_node(export_ctx, node.input("Detail")), name="detail")
    xml_node.add_child(export_node(export_ctx, node.input("Roughness")), name="roughness")
    xml_node.add_child(export_node(export_ctx, node.input("Lacunarity")), name="lacunarity")
    xml_node.add_child(export_node(export_ctx, node.input("Distortion")), name="distortion")

    return xml_node

def _export_separate_xyz_node(export_ctx: ExportContext, input: RMInput, exposure):
    factor = exposure or 1
    node = input.linked_node()
    output = input.linked_output()
    blnode = node.bl_node

    if output == "X":
        component="r"
    elif output == "Y":
        component="g"
    elif output == "Z":
        component="b"
    else:
        component=None
        export_ctx.error(f"[ShaderNodeSeparateXYZ] Connected to unrecognized output with name {output}")

    xml_node = XMLNode("texture", type="separate", component=component, convert="NONE")

    xml_node.add_child(export_node(export_ctx, node.input("Vector")))

    return xml_node

def _export_separate_color_node(export_ctx: ExportContext, input: RMInput, exposure):
    factor = exposure or 1
    node = input.linked_node()
    output = input.linked_output()
    blnode = node.bl_node

    convert_to = blnode.mode

    if convert_to == "RGB":
        convert_to = "NONE"

    if output == "Red":
        component="r"
    elif output == "Green":
        component="g"
    elif output == "Blue":
        component="b"
    else:
        component=None
        export_ctx.error(f"[ShaderNodeSeparateXYZ] Connected to unrecognized output with name {output}")

    xml_node = XMLNode("texture", type="separate", component=component, convert=convert_to)

    xml_node.add_child(export_node(export_ctx, node.input("Color")))

    return xml_node

def _export_ambient_occlusion(export_ctx: ExportContext, input: RMInput, exposure=None):
    node = input.linked_node().bl_node
    default_value = [ v for v in node.inputs[0].default_value[0:3] ]
    return XMLNode("texture", type="ao", distance=input.linked_node().bl_node.inputs[1].default_value, samples=input.linked_node().bl_node.samples, color=str_flat_array(default_value))

def _export_wave_texture(export_ctx: ExportContext, input: RMInput, exposure):
    factor = exposure or 1
    node = input.linked_node()
    output = input.linked_output()
    blnode = node.bl_node

    if blnode.wave_type == 'BANDS':
        xml_node = XMLNode("texture", type="wave", wavetype=blnode.wave_type, direction=blnode.bands_direction, waveprofile=blnode.wave_profile)
    else:
        export_ctx.error(f"Wave type '{blnode.wave_type}' is currently not supported.")
        return _export_fallback()

    xml_node.add_child(export_node(export_ctx, node.input("Scale")), name="scale")
    xml_node.add_child(export_node(export_ctx, node.input("Distortion")), name="distortion")
    xml_node.add_child(export_node(export_ctx, node.input("Detail")), name="detail")
    xml_node.add_child(export_node(export_ctx, node.input("Detail Scale")), name="detailscale")
    xml_node.add_child(export_node(export_ctx, node.input("Detail Roughness")), name="detailroughness")
    xml_node.add_child(export_node(export_ctx, node.input("Phase Offset")), name="phaseoffset")

    return xml_node

_node_handlers: dict[str, any] = {
    "ShaderNodeTexImage": _export_image_texture,
    "ShaderNodeTexEnvironment": _export_env_image_texture,
    "ShaderNodeValue": _export_scalar_value,
    "ShaderNodeRGB": _export_rgb_value,
    "ShaderNodeValToRGB": _export_val_to_rgb,
    "ShaderNodeNormalMap": _export_normal_map,
    "ShaderNodeParticleInfo": _export_particle_info,
    "ShaderNodeHairInfo": _export_curves_info,
    "ShaderNodeObjectInfo": _export_object_info,
    "ShaderNodeNewGeometry": _export_geometry_info,
    "ShaderNodeTexChecker": _export_checkerboard_texture,
    "ShaderNodeHueSaturation": _export_hue_saturation,
    "ShaderNodeTexNoise": _export_noise_texture,
    "ShaderNodeInvert": _export_invert,
    "ShaderNodeSeparateXYZ": _export_separate_xyz_node,
    "ShaderNodeSeparateColor": _export_separate_color_node,
    "ShaderNodeBrightContrast": _export_brightness_contrast,
    "ShaderNodeMath": _export_math,
    "ShaderNodeMapRange": _export_map_range,
    "ShaderNodeMix": _export_mix,
    "ShaderNodeFresnel": _export_fresnel,
    "ShaderNodeCombineXYZ": _export_combine_xyz,
    "ShaderNodeAmbientOcclusion": _export_ambient_occlusion,
    "ShaderNodeTexWave": _export_wave_texture,
}

def export_node(export_ctx: ExportContext, input: RMInput, exposure=None) -> XMLNode:
    if not input.is_linked():
        return _export_default(export_ctx, input, exposure)
    
    node = input.linked_node()
    for (typename, handler) in _node_handlers.items():
        if hasattr(bpy.types, typename) and isinstance(node.bl_node, getattr(bpy.types, typename)):
            return export_ctx.registry.export(Token(node.bl_node.as_pointer(), output=input.linked_output()),
                                   lambda: handler(export_ctx, input, exposure))
    
    export_ctx.error(f"Material {node.node_graph.name} has a node of type {type(node.bl_node).__name__} which is not supported")
    return _export_fallback()
