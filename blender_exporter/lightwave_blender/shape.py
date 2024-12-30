import os
import bmesh
import bpy
import traceback

from .addon_preferences import get_prefs
from .export_context import ExportContext
from .xml_node import XMLNode
from .materials import export_material
from .node_graph import RMInput, RMNode, RMNodeGraph
from .node import export_node
from .utils import LoggingBlock, logger


def get_shape_name_base(obj, inst):
    modifiers = [mod.type for mod in obj.original.modifiers]
    has_nodes = 'NODES' in modifiers

    if has_nodes:
        # Not sure how to ensure shapes with nodes are handled as uniques
        # TODO: We better join them by material
        id = hex(inst.random_id).replace("0x", "").replace('-', 'M').upper()
        return f"{obj.name}_{id}"

    try:
        return f"{obj.data.name}-shape"
    except:
        return f"{obj.original.data.name}-shape"  # We use the original mesh name!


def _shape_name_material(name, mat_id):
    return f"_m_{mat_id}_{name}"

def _export_bmesh_by_material(export_ctx: ExportContext, me, suffix: str) -> list[(str, str)]:
    mat_count = len(me.materials)
    shapes = []

    def _export_for_mat(mat_id, abs_filepath):
        from .mesh import mesh_save

        bm = bmesh.new()
        bm.from_mesh(me)

        # remove faces with other materials
        if mat_count > 1:
            for f in bm.faces:
                # Remove irrelevant faces
                # Special case: Assign invalid material indices to the last material 
                if f.material_index != mat_id and not ((f.material_index < 0 or f.material_index >= mat_count) and mat_id == mat_count-1):
                    bm.faces.remove(f)

        if len(bm.verts) == 0 or len(bm.faces) == 0 or not bm.is_valid:
            bm.free()
            return False

        # Make sure all faces are convex
        bmesh.ops.connect_verts_concave(bm, faces=bm.faces)
        bmesh.ops.triangulate(bm, faces=bm.faces)

        bm.normal_update()

        mesh_save(
            filepath=abs_filepath,
            bm=bm,
            )

        bm.free()
        return True
    
    if mat_count == 0:
        # special case if the mesh has no slots available
        mat_count = 1
    
    for mat_id in range(0, mat_count):
        with LoggingBlock(f"exporting mesh for material with id {mat_id}") as lb:
            shape_name = me.name if mat_count <= 1 else _shape_name_material(me.name, mat_id)
            rel_filepath = os.path.join(get_prefs().mesh_dir_name, shape_name + suffix + ".bob")
            abs_filepath = os.path.join(export_ctx.path, rel_filepath)

            if os.path.exists(abs_filepath) and not export_ctx.settings.overwrite_existing_meshes:
                # file is already exported
                pass
            elif _export_for_mat(mat_id, abs_filepath):
                # export successful
                pass
            else:
                # export failed
                continue
            
            shapes.append(rel_filepath.replace('\\', '/')) # Ensure the shape path is not using \ to keep the xml valid
    
    return shapes

def export_shape(depsgraph: bpy.types.Depsgraph, export_ctx: ExportContext, obj: bpy.types.Object, suffix: str) -> list[XMLNode]:
    obj_type = obj.type
    try:
        me = obj.to_mesh(preserve_all_data_layers=False, depsgraph=depsgraph)
    except RuntimeError as e:
        if obj_type not in {'MESH'}:
            export_ctx.warn(f"Dealing with an unsupported object of type {obj_type}")
            return [XMLNode("shape", type=f"{obj_type.lower()}")]
        # export_ctx.error(f"Could not convert to mesh: {str(e)}")
        # return []

    shapes = _export_bmesh_by_material(export_ctx, me, suffix)
    obj.to_mesh_clear()

    result = []
    for shape_index in range(len(shapes)):
        filepath = shapes[shape_index]
        
        if shape_index >= len(obj.material_slots):
            continue
        
        shape_mat = obj.material_slots[shape_index].material

        shape_is_displaced = False
        displacement = XMLNode("shape", type="displaced", filename=filepath)
        if shape_mat.use_nodes:
            try:
                node_graph = RMNodeGraph(export_ctx, shape_mat.name_full, shape_mat.node_tree)
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
                            input = rm_node.input("Displacement")
                            if input.is_linked():
                                displacement_node = input.linked_node()
                                print(displacement_node)
                                typename = "ShaderNodeDisplacement"
                                if hasattr(bpy.types, typename) and isinstance(displacement_node.bl_node, getattr(bpy.types, typename)):
                                    displacement.add_child(export_node(export_ctx, displacement_node.input("Height")), name="height")
                                    displacement.add_child(export_node(export_ctx, displacement_node.input("Midlevel")), name="offset")
                                    displacement.add_child(export_node(export_ctx, displacement_node.input("Scale")), name="scale")
                                    shape_is_displaced = True
            except Exception as e:
                print(f"failed to evaluate displacement for material {shape_mat.name}")
                print(e)
                traceback.print_exc()
        
        if (shape_is_displaced):
            result.append(displacement)
        else:
            result.append(XMLNode("shape", type="mesh", filename=filepath))

    return result
