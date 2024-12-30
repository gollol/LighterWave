import bpy
from .utils import *
from .defaults import *
from .export_context import ExportContext
from .xml_node import XMLNode
from .node_graph import RMNodeGraph, RMInput, RMNode
from .node import export_node
from .materials import _is_black
from .transform import export_transform_node


def export_world_background(export_ctx: ExportContext):
    with LoggingBlock("exporting world background") as lb:
        if not export_ctx.scene.world:
            return []

        # Export basic world if no shading nodes are given
        if not export_ctx.scene.world.node_tree:
            if export_ctx.scene.world.color[0] > 0 or export_ctx.scene.world.color[1] > 0 or export_ctx.scene.world.color[2] > 0:
                bg = XMLNode("light", type="simpleenvmap") #For constant environment maps this should always be better 
                bg.add("texture", type="constant",
                    value=str_flat_array(export_ctx.scene.world.color))
                return [bg]
            else:
                return []  # No world

        node_graph = RMNodeGraph(
            export_ctx, export_ctx.scene.world.name_full, export_ctx.scene.world.node_tree)
        node_graph.inline_node_groups_recursively()
        node_graph.remove_reroute_nodes()
        node_graph.remove_muted_nodes()
        node_graph.remove_layout_nodes()

        for (node_name, rm_node) in node_graph.nodes.items():
            node = rm_node.bl_node
            if isinstance(node, bpy.types.ShaderNodeOutputWorld):
                if node.is_active_output:
                    return _export_world(export_ctx, rm_node.input("Surface"))

        return []  # No active output

def _get_bg_xml_node(type: str, sample_type: str) -> XMLNode:
    if type == "SIMPLE":
        bg = XMLNode("light", type="simpleenvmap")
    elif type == "WARP":
        bg = XMLNode("light", type="warp2devnmap", )
        if sample_type == "MIS":
            compensation = 0.5
        else:
            compensation = 0.0
        bg.add("float", name="mis_compensation_strength", value=compensation)
    else:
        bg = XMLNode("light", type="envmap")
    return bg

def _export_background(export_ctx: ExportContext, bsdf_node: RMNode):
    color = bsdf_node.input("Color")
    strength = bsdf_node.input("Strength")

    if not color.is_linked():
        if not color.has_value() or _is_black(color.value):
            return []

    emission_scale = 1
    if strength.is_linked():
        export_ctx.error(
            "Only constant values for emission strength are supported")
    elif strength.has_value():
        emission_scale = float(strength.value)
    if emission_scale == 0:
        return []

    bg = _get_bg_xml_node(export_ctx.settings.environment_map_sampling_mode, export_ctx.settings.sampling_mode)

    bg.add_child(export_node(
        export_ctx, color, exposure=emission_scale))

    transforms = []
    transforms += [XMLNode("matrix", value=str_flat_matrix(ENVIRONMENT_MAP_TRANSFORM))]
    transforms += export_transform_node(export_ctx, color)

    bg.add("transform").add_children(transforms)
    return [bg]


_world_handlers: dict[str, any] = {
    "ShaderNodeEmission": _export_background,
    "ShaderNodeBackground": _export_background,
}


def _export_world(export_ctx: ExportContext, input: RMInput):
    if not input.is_linked():
        return []  # black world

    world_node = input.linked_node()

    if world_node is None:
        export_ctx.error(f"Material {input.node_graph.name} has no valid node")
        return []

    result = []
    for (typename, handler) in _world_handlers.items():
        if hasattr(bpy.types, typename) and isinstance(world_node.bl_node, getattr(bpy.types, typename)):
            result += handler(export_ctx, world_node)
            break
    else:
        # treat as background
        transforms = []
        transforms += [XMLNode("matrix", value=str_flat_matrix(ENVIRONMENT_MAP_TRANSFORM))]
        transforms += export_transform_node(export_ctx, input)

        bg = _get_bg_xml_node(export_ctx.settings.environment_map_sampling_mode, export_ctx.settings.sampling_mode)
        bg.add_child(export_node(export_ctx, input))
        bg.add("transform").add_children(transforms)
        result.append(bg)

    return result
