import bpy
from .utils import *
from .defaults import *
from .export_context import ExportContext
from .xml_node import XMLNode
from .node_graph import RMInput, RMNode
from typing import Callable

import math

def create_motion_blur_transform(export_ctx: ExportContext, instance: bpy.types.DepsgraphObjectInstance, transform_fn: Callable[[bpy.types.DepsgraphObjectInstance],XMLNode]) -> XMLNode:
    """Creates an animated transform for motion blur using the stub renderer"""

    #NOTE: This is not the best way to do it, as it can be slow,
    #      but was easy to integrate into architecture
    
    ani_trans_node = XMLNode("transform", type="animated")
    ani_trans_node.add_child(transform_fn(instance))

    has_animation = (instance.object is not None) and (instance.object.animation_data is not None)
    
    # Handle motion blur
    if has_animation:
        if (export_ctx.engine is not None) and export_ctx.scene.render.use_motion_blur:
            original_frame = export_ctx.scene.frame_current
            original_subframe = export_ctx.scene.frame_subframe
            export_ctx.engine.frame_set(original_frame, export_ctx.scene.render.motion_blur_shutter)

            # Get transform for next step for motion blur
            ani_trans_node.add_child(transform_fn(instance))
            
            # Reset frame
            export_ctx.engine.frame_set(original_frame, original_subframe)
        
        return ani_trans_node
    else:
        return transform_fn(instance)

def _extract_vector(export_ctx: ExportContext, input: RMInput):
    if input.is_linked():
        export_ctx.error(
            "Only constant values for transformations are supported")

    if not input.has_value():
        return [0, 0, 0]

    if isinstance(input.value, float):
        return [input.value, input.value, input.value]
    else:
        return input.value


def _export_sca(v):
    if any((c != 1 for c in v[0:3])):
        return [XMLNode("scale", value=str_flat_array([ 1/x for x in v[0:3] ]))]
    else:
        return []


def _export_rot(v):
    res = []
    if v[2] != 0:
        res += [XMLNode("rotate", axis="0,0,1", angle=-v[2])]
    if v[1] != 0:
        res += [XMLNode("rotate", axis="0,1,0", angle=-v[1])]
    if v[0] != 0:
        res += [XMLNode("rotate", axis="1,0,0", angle=-v[0])]
    return res


def _export_tra(v):
    if any((c != 0 for c in v[0:3])):
        return [XMLNode("translate", value=str_flat_array([ -x for x in v[0:3] ]))]
    else:
        return []


def _export_vector_mapping(export_ctx: ExportContext, node: RMNode):
    sca = _extract_vector(export_ctx, node.input("Scale"))
    rot = _extract_vector(
        export_ctx, node.input("Rotation"))
    loc = _extract_vector(
        export_ctx, node.input("Location"))

    rot = [r * 180 / math.pi for r in rot]  # rad to deg

    transforms = export_transform_node(export_ctx, node.input("Vector"))

    if node.bl_node.vector_type == 'POINT':
        transforms += _export_tra(loc)
        transforms += _export_rot(rot)
        transforms += _export_sca(sca)
    elif node.bl_node.vector_type == 'TEXTURE':
        transforms += _export_sca([1.0/v for v in sca])
        transforms += list(reversed(_export_rot([-v for v in rot])))
        transforms += _export_tra([-v for v in loc])
    elif node.bl_node.vector_type == 'NORMAL':
        # No idea why like this
        transforms += _export_rot(rot)
        transforms += _export_sca([1.0/v for v in sca])
        # Usually the output has to be normalized here
    elif node.bl_node.vector_type == 'VECTOR':
        transforms += _export_rot(rot)
        transforms += _export_sca(sca)
    else:
        export_ctx.error(
            f"Material {node.node_graph.name} has a mapping of type {node.bl_node.vector_type} which is not supported")
    return transforms


def _export_vector_rotate(export_ctx: ExportContext, node: RMNode):
    center = _extract_vector(export_ctx, node.input("Center"))

    transforms = export_transform_node(export_ctx, node.input("Vector"))
    transforms += _export_tra([-c for c in center])
    if node.bl_node.rotation_type == "EULER_XYZ":
        rot = _extract_vector(export_ctx, node.input("Rotation"))

        rot = [r * 180 / math.pi for r in rot]  # rad to deg

        if node.bl_node.invert:
            transforms += list(reversed(_export_rot([-a for a in rot])))
        else:
            transforms += _export_rot(rot)
    else:
        if node.input("Angle").has_value():
            angle = node.input("Angle").value * 180 / math.pi  # rad to deg
        else:
            return []

        if node.bl_node.invert:
            angle = -angle

        if angle == 0:
            return []

        if node.bl_node.rotation_type == "AXIS_ANGLE":
            axis = _extract_vector(export_ctx, node.input("Axis"))
            transforms += [XMLNode("rotate",
                                   axis=str_flat_array(axis), angle=-angle)]
        elif node.bl_node.rotation_type == "X_AXIS":
            transforms += [XMLNode("rotate", axis="1,0,0", angle=-angle)]
        elif node.bl_node.rotation_type == "Y_AXIS":
            transforms += [XMLNode("rotate", axis="0,1,0", angle=-angle)]
        elif node.bl_node.rotation_type == "Z_AXIS":
            transforms += [XMLNode("rotate", axis="0,0,1", angle=-angle)]

    transforms += _export_tra(center)
    return transforms


def _export_vector_forward(export_ctx: ExportContext, node: RMNode):
    if node.input("Vector").is_linked():
        return export_transform_node(export_ctx, node.input("Vector"))
    else:
        return []


def _export_tex_coord(export_ctx: ExportContext, node: RMNode):
    return []  # Ignore for now


_node_handlers: dict[str, any] = {
    "ShaderNodeTexCoord": _export_tex_coord,
    "ShaderNodeMapping": _export_vector_mapping,
    "ShaderNodeVectorRotate": _export_vector_rotate,
    "ShaderNodeTexImage": _export_vector_forward,
    "ShaderNodeTexEnvironment": _export_vector_forward,
}


def export_transform_node(export_ctx: ExportContext, input: RMInput) -> XMLNode:
    if not input.is_linked():
        return []

    node = input.linked_node()

    result = []
    for (typename, handler) in _node_handlers.items():
        if hasattr(bpy.types, typename) and isinstance(node.bl_node, getattr(bpy.types, typename)):
            result += handler(export_ctx, node)
            break
    else:
        export_ctx.error(
            f"Node graph {node.node_graph.name} has a node of type {type(node.bl_node).__name__} which is not supported")
    return result
