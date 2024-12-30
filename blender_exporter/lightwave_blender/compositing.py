import bpy

from .xml_node import XMLNode
from .node_graph import RMInput, RMNodeGraph
from .export_context import ExportContext
from .registry import Token
from .utils import *

#
#   Main export functions
#

def export_render_passes(export_ctx: ExportContext, view_layer: bpy.types.ViewLayer, name_prefix: str) -> list[XMLNode]:
    """Export integrator nodes for all the different selected render passes"""

    if getattr(export_ctx.scene, 'cycles', None) is not None and bpy.context.engine == 'CYCLES':
        cycles = export_ctx.scene.cycles
        max_depth = cycles.max_bounces + 1
        clamp = max(cycles.sample_clamp_direct, cycles.sample_clamp_indirect)
        spp = cycles.samples
        use_denoising = cycles.use_denoising
        denoising_store_passes = view_layer.cycles.denoising_store_passes
    else:
        max_depth = 10
        clamp = 0
        spp = 64
        use_denoising = False
        denoising_store_passes = False

    def _img_name(pass_name: str):
        return get_image_name(view_layer, pass_name, name_prefix)

    result = []

    pt = XMLNode("integrator", type="pathtracer", depth=max_depth, samplingmode=export_ctx.settings.sampling_mode)
    pt.add("ref", id="scene")
    pt_img = XMLNode("image", id=_img_name("pathtracer"))
    pt.add_child(export_ctx.registry.export(Token("integrator", id=_img_name("pathtracer")), lambda: pt_img))
    pt.add("sampler", type=export_ctx.settings.sampler_type.lower(), count=spp)
    result.append(pt)

    if view_layer.use_pass_normal or denoising_store_passes:
        norm = XMLNode("integrator", type="aov", variable="normals")
        norm.add("ref", id="scene")
        norm_img = XMLNode("image", id=_img_name("normals"))
        norm.add_child(export_ctx.registry.export(Token("integrator", id=_img_name("normals")), lambda: norm_img))
        norm.add("sampler", type=export_ctx.settings.sampler_type.lower(), count=spp)
        result.append(norm)
    
    if denoising_store_passes:
        alb = XMLNode("integrator", type="aov", variable="albedo")
        alb.add("ref", id="scene")
        alb_img = XMLNode("image", id=_img_name("albedo"))
        alb.add_child(export_ctx.registry.export(Token("integrator", id=_img_name("albedo")), lambda: alb_img))
        alb.add("sampler", type=export_ctx.settings.sampler_type.lower(), count=spp)
        result.append(alb)
    
    if view_layer.use_pass_z:
        z = XMLNode("integrator", type="aov", variable="distance")
        z.add("ref", id="scene")
        z_img = XMLNode("image", id=_img_name("distance"))
        z.add_child(export_ctx.registry.export(Token("integrator", id=_img_name("distance")), lambda: z_img))
        z.add("sampler", type=export_ctx.settings.sampler_type.lower(), count=spp)
        result.append(z)
    
    if view_layer.use_pass_uv:
        uv = XMLNode("integrator", type="aov", variable="uv")
        uv.add("ref", id="scene")
        uv_img = XMLNode("image", id=_img_name("uv"))
        uv.add_child(export_ctx.registry.export(Token("integrator", id=_img_name("uv")), lambda: uv_img))
        uv.add("sampler", type=export_ctx.settings.sampler_type.lower(), count=spp)
        result.append(uv)
    
    if view_layer.use_pass_ambient_occlusion:
        ao = XMLNode("integrator", type="ao", distance=1)
        ao.add("ref", id="scene")
        ao_img = XMLNode("image", id=_img_name("ao"))
        ao.add_child(export_ctx.registry.export(Token("integrator", id=_img_name("ao")), lambda: ao_img))
        ao.add("sampler", type=export_ctx.settings.sampler_type.lower(), count=spp)
        result.append(ao)
    
    return result

def export_composition(export_ctx: ExportContext) -> list[XMLNode]:
    """Export composition graph as postprocess nodes in Lightwave"""

    with LoggingBlock(f"exporting composition") as lb:
        scene = export_ctx.scene
        name_prefix = f"{export_ctx.settings.scene_name_prefix.strip()}_" if len(export_ctx.settings.scene_name_prefix.strip()) > 0 else ""

        result: list[XMLNode] = []

        # Get render passes for each view layer, as they have to be rendered separately
        for view_layer in export_ctx.scene.view_layers:
            if not view_layer.use or not export_ctx.is_view_layer_selected(view_layer):
                lb.debug(f"Skipping view layer '{view_layer.name}' as it is disabled or not selected.")
                continue

            result = result + export_render_passes(export_ctx, view_layer, name_prefix)

        if not scene.use_nodes:
            return result

        # Cleanup node graph from unnecessary constructs
        node_graph = RMNodeGraph(export_ctx, "Composition", scene.node_tree)
        node_graph.inline_node_groups_recursively()
        node_graph.remove_reroute_nodes()
        node_graph.remove_muted_nodes()
        node_graph.remove_layout_nodes()
        
        for (node_name, rm_node) in node_graph.nodes.items():
            node = rm_node.bl_node

            if isinstance(node, bpy.types.CompositorNodeComposite):
                img = export_compositor_node(export_ctx, result, name_prefix, rm_node.input("Image"))

                #TODO: Multiple Composite Nodes??

        return result

#
#   Node exporting
#

def _export_fallback_image() -> XMLNode:
    """Export a fallback image of constant black color"""
    return XMLNode("image", color="0")

def _export_default_image(export_ctx: ExportContext, input: RMInput) -> XMLNode:
    """Export default image with constant color, if present"""

    if not input.has_value():
        return _export_fallback_image()

    value = input.value
    if isinstance(value, list):
        value = value[0:3]
    
    return XMLNode("image", color=str_flat_array(value))

_render_pass_label_id_map = {
    "Image":            "pathtracer",
    "Depth":            "distance",
    "UV":               "uv",
    "AO":               "ao",
    "Normal":           "normals",
    "Denoising Normal": "normals",
    "Denoising Albedo": "albedo",
    "Denoising Depth":  "distance",
}

def _export_render_layer_compositor_node(export_ctx: ExportContext, result: list[XMLNode], name_prefix: str, input: RMInput) -> XMLNode:
    node = input.linked_node()
    output = input.linked_output()
    bl_node = node.bl_node

    if not (output in _render_pass_label_id_map):
        export_ctx.error(f"Output '{input.linked_output()}' of 'Render Layers' node is not supported!")
        return _export_fallback_image()
    
    view_layer = export_ctx.scene.view_layers[bl_node.layer]

    if not export_ctx.is_view_layer_selected(view_layer):
        export_ctx.error(f"Render Layer Composition node used for layer '{bl_node.layer}', which is not selected for exporting!")
        return _export_fallback_image()
    
    return export_ctx.registry.export(Token("integrator", id=get_image_name(view_layer, _render_pass_label_id_map[output], name_prefix)), _export_fallback_image)

def _export_exposure_compositor_node(export_ctx: ExportContext, result: list[XMLNode], name_prefix: str, input: RMInput) -> XMLNode:
    node = input.linked_node()
    output = input.linked_output()
    bl_node = node.bl_node

    input_img = export_compositor_node(export_ctx, result, name_prefix, node.input("Image"))
    input_exp = export_compositor_node(export_ctx, result, name_prefix, node.input("Exposure"))

    xml_node = XMLNode("postprocess", type="exposure")
    xml_node.add_child(input_img, name="input")
    xml_node.add_child(input_exp, name="exposure")

    output_img = XMLNode("image", id=export_ctx.registry.make_unique_id(name_prefix + "compositor_exposure"))
    xml_node.add_child(export_ctx.registry.export(Token("composition_node", p=bl_node.as_pointer(), o=input.linked_output()), lambda: output_img), name="output")

    result.append(xml_node)
    return output_img

def _export_denoise_compositor_node(export_ctx: ExportContext, result: list[XMLNode], name_prefix: str, input: RMInput) -> XMLNode:
    node = input.linked_node()
    output = input.linked_output()
    bl_node = node.bl_node

    input_img = export_compositor_node(export_ctx, result, name_prefix, node.input("Image"))
    input_norm = export_compositor_node(export_ctx, result, name_prefix, node.input("Normal"))
    input_alb = export_compositor_node(export_ctx, result, name_prefix, node.input("Albedo"))

    xml_node = XMLNode("postprocess", type="denoise", hdr=bl_node.use_hdr, prefilter=bl_node.prefilter)
    xml_node.add_child(input_img, name="input")
    xml_node.add_child(input_norm, name="normals")
    xml_node.add_child(input_alb, name="albedo")

    output_img = XMLNode("image", id=export_ctx.registry.make_unique_id(name_prefix + "compositor_denoise"))
    xml_node.add_child(export_ctx.registry.export(Token("composition_node", p=bl_node.as_pointer(), o=input.linked_output()), lambda: output_img), name="output")

    result.append(xml_node)
    return output_img

def _export_tonemap_compositor_node(export_ctx: ExportContext, result: list[XMLNode], name_prefix: str, input: RMInput) -> XMLNode:
    node = input.linked_node()
    output = input.linked_output()
    bl_node = node.bl_node

    if (bl_node.tonemap_type != 'RH_SIMPLE'):
        export_ctx.error(f"'{bl_node.name}' compositor node uses a tonemapping type other than 'Rh Simple', which is not supported!")
        return _export_default_image(export_ctx, input)

    input_img = export_compositor_node(export_ctx, result, name_prefix, node.input("Image"))

    # Currently, only "ReinhardSimple" is supported, so this is hardcoded
    xml_node = XMLNode("postprocess", type="tonemapping", key=bl_node.key, tmmode="ReinhardSimple")
    xml_node.add_child(input_img, name="input")

    output_img = XMLNode("image", id=export_ctx.registry.make_unique_id(name_prefix + "compositor_tonemapping"))
    xml_node.add_child(export_ctx.registry.export(Token("composition_node", p=bl_node.as_pointer(), o=input.linked_output()), lambda: output_img), name="output")

    result.append(xml_node)
    return output_img

def _export_bright_contrast_compositor_node(export_ctx: ExportContext, result: list[XMLNode], name_prefix: str, input: RMInput) -> XMLNode:
    node = input.linked_node()
    output = input.linked_output()
    bl_node = node.bl_node

    input_img = export_compositor_node(export_ctx, result, name_prefix, node.input("Image"))
    input_bright = export_compositor_node(export_ctx, result, name_prefix, node.input("Bright"))
    input_contrast = export_compositor_node(export_ctx, result, name_prefix, node.input("Contrast"))

    xml_node = XMLNode("postprocess", type="brightness_contrast")
    xml_node.add_child(input_img, name="input")
    xml_node.add_child(input_bright, name="brightness")
    xml_node.add_child(input_contrast, name="contrast")

    output_img = XMLNode("image", id=export_ctx.registry.make_unique_id(name_prefix + "compositor_brightness_contrast"))
    xml_node.add_child(export_ctx.registry.export(Token("composition_node", p=bl_node.as_pointer(), o=input.linked_output()), lambda: output_img), name="output")

    result.append(xml_node)
    return output_img

def _export_val_to_rgb_compositor_node(export_ctx: ExportContext, result: list[XMLNode], name_prefix: str, input: RMInput) -> XMLNode:
    node = input.linked_node()
    output = input.linked_output()
    bl_node = node.bl_node

    color_ramp = bl_node.color_ramp

    input_img = export_compositor_node(export_ctx, result, name_prefix, node.input("Image"))
    input_fac = export_compositor_node(export_ctx, result, name_prefix, node.input("Fac"))

    # Interpolation type is different for different color modes, so only export the used one
    if color_ramp.color_mode == 'RGB':
        xml_node = XMLNode("postprocess", type="colorramp", colorMode=color_ramp.color_mode, interpolationType=color_ramp.interpolation)
    else:
        xml_node = XMLNode("postprocess", type="colorramp", colorMode=color_ramp.color_mode, hueInterpolation=color_ramp.hue_interpolation)

    xml_node.add_child(input_img, name="input")
    xml_node.add_child(input_fac, name="factor")

    for element in color_ramp.elements.values():
        element_xml_node = XMLNode("colorRampElement", alpha=element.alpha, position=element.position)
        element_xml_node.add_child(XMLNode("color", value=str_flat_array(list(element.color[:3]))), name="color")

        xml_node.add_child(element_xml_node)

    output_img = XMLNode("image", id=export_ctx.registry.make_unique_id(name_prefix + "compositor_colorramp"))
    xml_node.add_child(export_ctx.registry.export(Token("composition_node", p=bl_node.as_pointer(), o=input.linked_output()), lambda: output_img), name="output")

    result.append(xml_node)
    return output_img

def _export_hue_sat_compositor_node(export_ctx: ExportContext, result: list[XMLNode], name_prefix: str, input: RMInput) -> XMLNode:
    node = input.linked_node()
    output = input.linked_output()
    bl_node = node.bl_node

    input_img = export_compositor_node(export_ctx, result, name_prefix, node.input("Image"))
    input_hue = export_compositor_node(export_ctx, result, name_prefix, node.input("Hue"))
    input_sat = export_compositor_node(export_ctx, result, name_prefix, node.input("Saturation"))
    input_val = export_compositor_node(export_ctx, result, name_prefix, node.input("Value"))
    #input_fac = export_compositor_node(export_ctx, result, name_prefix, node.input("Fac"))

    xml_node = XMLNode("postprocess", type="hue_saturation_value")
    xml_node.add_child(input_img, name="input")
    xml_node.add_child(input_hue, name="hue")
    xml_node.add_child(input_sat, name="saturation")
    xml_node.add_child(input_val, name="value")
    #xml_node.add_child(input_fac, name="factor")

    output_img = XMLNode("image", id=export_ctx.registry.make_unique_id(name_prefix + "compositor_hue_saturation_value"))
    xml_node.add_child(export_ctx.registry.export(Token("composition_node", p=bl_node.as_pointer(), o=input.linked_output()), lambda: output_img), name="output")

    result.append(xml_node)
    return output_img

def _export_math_compositor_node(export_ctx: ExportContext, result: list[XMLNode], name_prefix: str, input: RMInput) -> XMLNode:
    node = input.linked_node()
    output = input.linked_output()
    bl_node = node.bl_node

    input_val1 = export_compositor_node(export_ctx, result, name_prefix, node.input("Value"))
    input_val2 = export_compositor_node(export_ctx, result, name_prefix, node.input("Value_001"))
    input_val3 = export_compositor_node(export_ctx, result, name_prefix, node.input("Value_002"))

    xml_node = XMLNode("postprocess", type="math", clamping=bl_node.use_clamp, operation=bl_node.operation.lower())
    xml_node.add_child(input_val1, name="input")
    xml_node.add_child(input_val2, name="b_image")
    #xml_node.add_child(input_val3, name="c_image")

    output_img = XMLNode("image", id=export_ctx.registry.make_unique_id(name_prefix + "compositor_math"))
    xml_node.add_child(export_ctx.registry.export(Token("composition_node", p=bl_node.as_pointer(), o=input.linked_output()), lambda: output_img), name="output")

    result.append(xml_node)
    return output_img

def _export_mix_rgb_compositor_node(export_ctx: ExportContext, result: list[XMLNode], name_prefix: str, input: RMInput) -> XMLNode:
    node = input.linked_node()
    output = input.linked_output()
    bl_node = node.bl_node

    input_img1 = export_compositor_node(export_ctx, result, name_prefix, node.input("Image"))
    input_img2 = export_compositor_node(export_ctx, result, name_prefix, node.input("Image_001"))
    input_fac = export_compositor_node(export_ctx, result, name_prefix, node.input("Fac"))

    xml_node = XMLNode("postprocess", type="mixColor", clamping=bl_node.use_clamp, operation=bl_node.blend_type.lower())
    xml_node.add_child(input_img1, name="input")
    xml_node.add_child(input_img2, name="b_image")
    xml_node.add_child(input_fac, name="factor")

    output_img = XMLNode("image", id=export_ctx.registry.make_unique_id(name_prefix + "compositor_mixColor"))
    xml_node.add_child(export_ctx.registry.export(Token("composition_node", p=bl_node.as_pointer(), o=input.linked_output()), lambda: output_img), name="output")

    result.append(xml_node)
    return output_img

_compositor_node_handlers: dict[str, any] = {
    "CompositorNodeRLayers": _export_render_layer_compositor_node,
    "CompositorNodeExposure": _export_exposure_compositor_node,
    "CompositorNodeDenoise": _export_denoise_compositor_node,
    "CompositorNodeTonemap": _export_tonemap_compositor_node,
    "CompositorNodeBrightContrast": _export_bright_contrast_compositor_node,
    "CompositorNodeValToRGB": _export_val_to_rgb_compositor_node,
    "CompositorNodeHueSat": _export_hue_sat_compositor_node,
    "CompositorNodeMath": _export_math_compositor_node,
    "CompositorNodeMixRGB": _export_mix_rgb_compositor_node,
}

def export_compositor_node(export_ctx: ExportContext, result: list[XMLNode], name_prefix: str, input: RMInput) -> XMLNode:
    """Export a single compositor node based on node type and using the scene registry for deduplication when possible"""
    if not input.is_linked():
        return _export_default_image(export_ctx, input)
    
    node = input.linked_node()

    for (typename, handler) in _compositor_node_handlers.items():
        if hasattr(bpy.types, typename) and isinstance(node.bl_node, getattr(bpy.types, typename)):
            xml_node = handler(export_ctx, result, name_prefix, input)
            return export_ctx.registry.export(Token("composition_node", p=node.bl_node.as_pointer(), o=input.linked_output()), lambda: xml_node)
    
    export_ctx.error(f"Composition has a node of type {type(node.bl_node).__name__} which is not supported")
    return _export_fallback_image()