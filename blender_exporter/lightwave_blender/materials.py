import bpy

from .export_context import ExportContext
from .node import export_node
from .utils import *
from .xml_node import XMLNode
from .node_graph import RMInput, RMNode, RMNodeGraph


def export_default_bsdf():
    node = XMLNode("bsdf", type="diffuse")
    node.add("texture", name="albedo", type="constant", value=0.8)
    return [node]

def _export_diffuse_bsdf(export_ctx: ExportContext, bsdf_node: RMNode):
    node = XMLNode("bsdf", type="diffuse")
    node.add_child(export_node(export_ctx, bsdf_node.input("Color")), name="albedo")
    return [node]

def _export_translucent_bsdf(export_ctx: ExportContext, bsdf_node: RMNode):
    node = XMLNode("bsdf", type="translucent")
    node.add_child(export_node(export_ctx, bsdf_node.input("Color")), name="albedo")
    return [node]

def _export_glass_bsdf_helper(export_ctx: ExportContext, bsdf_node: RMNode, has_reflectance: bool):
    has_roughness = bsdf_node.bl_node.distribution != 'SHARP'
    node = XMLNode("bsdf", type="roughdielectric" if has_roughness else "dielectric")
    node.add_child(export_node(export_ctx, bsdf_node.input("IOR")), name="ior")
    node.add_child(export_node(export_ctx, bsdf_node.input("Color")), name="transmittance")
    
    if has_reflectance:
        node.add_child(export_node(export_ctx, bsdf_node.input("Color")), name="reflectance")
    else:
        node.add("texture", name="reflectance", type="constant", value=0)
    
    if has_roughness:
        node.add_child(export_node(export_ctx, bsdf_node.input("Roughness")), name="roughness")
    return [node]

def _export_glass_bsdf(export_ctx: ExportContext, bsdf_node: RMNode):
    return _export_glass_bsdf_helper(export_ctx, bsdf_node, has_reflectance=True)

def _export_refraction_bsdf(export_ctx: ExportContext, bsdf_node: RMNode):
    return _export_glass_bsdf_helper(export_ctx, bsdf_node, has_reflectance=False)

def _export_transparent_bsdf(export_ctx: ExportContext, bsdf_node: RMNode):
    node = XMLNode("bsdf", type="dielectric")
    node.add("texture", name="ior", type="constant", value=1)
    node.add("texture", name="reflectance", type="constant", value=0)
    node.add_child(export_node(export_ctx, bsdf_node.input("Color")), name="transmittance")
    return [node]

def _export_glossy_bsdf(export_ctx: ExportContext, bsdf_node: RMNode):
    has_roughness = bsdf_node.input("Roughness").value > 0
    node = XMLNode("bsdf", type="roughconductor" if has_roughness else "conductor")
    node.add_child(export_node(export_ctx, bsdf_node.input("Color")), name="reflectance")
    if has_roughness:
        node.add_child(export_node(export_ctx, bsdf_node.input("Roughness")), name="roughness")
    return [node]

def _is_black(v: list):
    if isinstance(v, float):
        return v == 0
    return all((c == 0 for c in v[0:3]))

def _export_emission_helper(export_ctx: ExportContext, color: RMInput, strength: RMInput):
    if not color.is_linked():
        if not color.has_value() or _is_black(color.value):
            return []
    
    emission_scale = 1
    if strength.is_linked():
        export_ctx.error("Only constant values for emission strength are supported")
    elif strength.has_value():
        emission_scale = float(strength.value)
    if emission_scale == 0:
        return []

    emission = XMLNode("emission", type="lambertian")
    emission.add_child(export_node(export_ctx, color, exposure=emission_scale), name="emission")
    return [emission]

def _export_principled_bsdf(export_ctx: ExportContext, bsdf_node: RMNode):
    node = XMLNode("bsdf", type="principled")

    if (bsdf_node.bl_node.distribution != "GGX"):
        export_ctx.warn("Principled BSDF is using a different microfacet distribution than GGX.\nExported scene might look different!")

    for (lw_name, bl_name) in [
        ("baseColor", "Base Color"),
        ("metallic", "Metallic"),
        ("roughness", "Roughness"),
        ("ior", "IOR"),
        ("specularLevel", "Specular IOR Level"),
        ("specularTint", "Specular Tint"),
        ("transmission", "Transmission Weight"),
    ]:
        node.add_child(export_node(export_ctx, bsdf_node.input(bl_name)), name=lw_name)
    
    emission = _export_emission_helper(export_ctx, bsdf_node.input("Emission"), bsdf_node.input("Emission Strength"))
    return [node] + emission

def _export_emission(export_ctx: ExportContext, bsdf_node: RMNode):
    return _export_emission_helper(export_ctx, bsdf_node.input("Color"), bsdf_node.input("Strength"))

def _export_hair_bsdf(export_ctx: ExportContext, bsdf_node: RMNode):
    comp = bsdf_node.bl_node.component
    node = XMLNode("bsdf", type="hair", component=comp.lower())
    for (lw_name, bl_name) in [
        ("color", "Color"),
        ("offset", "Offset"),
        ("roughnessU", "RoughnessU"),
        ("roughnessV", "RoughnessV"),
        ("tangent", "Tangent"),
    ]:
        node.add_child(export_node(export_ctx, bsdf_node.input(bl_name)), name=lw_name)
    return [node]

def _export_mix_shader(export_ctx: ExportContext, bsdf_node: RMNode):
    node = XMLNode("bsdf", type="mix_shader")

    node.add_child(export_node(export_ctx, bsdf_node.input("Fac")), name="factor")

    node.add_children(_export_bsdf(export_ctx, bsdf_node.input("Shader")))
    node.add_children(_export_bsdf(export_ctx, bsdf_node.input("Shader_001")))

    return [node]

def _export_add_shader(export_ctx: ExportContext, bsdf_node: RMNode):
    node = XMLNode("bsdf", type="add_shader")

    node.add_children(_export_bsdf(export_ctx, bsdf_node.input("Shader")))
    node.add_children(_export_bsdf(export_ctx, bsdf_node.input("Shader_001")))

    return [node]


_bsdf_handlers: dict[str, any] = {
    "ShaderNodeBsdfDiffuse": _export_diffuse_bsdf,
    "ShaderNodeBsdfTranslucent": _export_translucent_bsdf,
    "ShaderNodeBsdfGlass": _export_glass_bsdf,
    "ShaderNodeBsdfRefraction": _export_refraction_bsdf,
    "ShaderNodeBsdfTransparent": _export_transparent_bsdf,
    "ShaderNodeBsdfAnisotropic": _export_glossy_bsdf,
    "ShaderNodeBsdfPrincipled": _export_principled_bsdf,
    "ShaderNodeBsdfHair": _export_hair_bsdf,
    "ShaderNodeEmission": _export_emission,
    "ShaderNodeBackground": _export_emission,
    "ShaderNodeMixShader": _export_mix_shader,
    "ShaderNodeAddShader": _export_add_shader
}

# @todo material type should be 'Material | World | Light'
def export_material(export_ctx: ExportContext, material: bpy.types.Material):
    with LoggingBlock(f"exporting material '{material.name}'") as lb:
        if not material.use_nodes:
            if isinstance(material, bpy.types.Light):
                return [] # TODO
            elif isinstance(material, bpy.types.World):
                return [] # TODO
            elif isinstance(material, bpy.types.Material):
                return export_default_bsdf() # TODO
            else:
                export_ctx.error("Unsupported use of node trees")
                return []
        
        try:
            node_graph = RMNodeGraph(export_ctx, material.name_full, material.node_tree)
            node_graph.inline_node_groups_recursively()
            node_graph.remove_reroute_nodes()
            node_graph.remove_muted_nodes()
            node_graph.remove_layout_nodes()

            for (node_name, rm_node) in node_graph.nodes.items():
                node = rm_node.bl_node
                if isinstance(node, bpy.types.ShaderNodeOutputMaterial) \
                or isinstance(node, bpy.types.ShaderNodeOutputWorld) \
                or isinstance(node, bpy.types.ShaderNodeOutputLight):
                    if node.is_active_output:
                        return _export_bsdf(export_ctx, rm_node.input("Surface"))
        except Exception as e:
            lb.error(f"Failed to export material {material.name}", exc_info=True)
        
        return [] # No active output

def _export_bsdf(export_ctx: ExportContext, input: RMInput):
    if not input.is_linked():
        return [] # black bsdf

    bsdf_node = input.linked_node()
    
    if bsdf_node is None:
        export_ctx.error(f"Material {input.node_graph.name} has no valid bsdf")
        return []
    
    result = []
    for (typename, handler) in _bsdf_handlers.items():
        if hasattr(bpy.types, typename) and isinstance(bsdf_node.bl_node, getattr(bpy.types, typename)):
            result += handler(export_ctx, bsdf_node)
            break
    else:
        # treat as emission
        emission = XMLNode("emission", type="lambertian")
        emission.add_child(export_node(export_ctx, input), name="emission")
        result.append(emission)

    if (normal := bsdf_node.input("Normal")).is_linked():
        normal_node = export_node(export_ctx, normal)
        normal_node.set_name("normal")
        result.append(normal_node)

    if (alpha := bsdf_node.input("Alpha")).is_linked() or (alpha.has_value() and alpha.value != 1):
        alpha_node = export_node(export_ctx, alpha)
        alpha_node.set_name("alpha")
        result.append(alpha_node)

    return result
