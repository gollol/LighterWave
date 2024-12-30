import bpy
import os

from .light import export_light
from .shape import export_shape
from .camera import stub_export_camera
from .materials import export_material, export_default_bsdf
from .utils import *
from .defaults import *
from .addon_preferences import get_prefs
from .xml_node import XMLNode, XMLRootNode
from .registry import SceneRegistry
from .export_context import ExportContext
from .world import export_world_background
from .fur import export_fur
from .compositing import export_composition
from datetime import datetime
from .transform import create_motion_blur_transform

class LightwaveStubRenderEngine(bpy.types.RenderEngine):
    """Stub rendering engine just used to get the render depsgraph for a view layer and export that"""

    # Set identifers and settings for stub renderer (as few capabilities as possible, since we just need the evaluated depsgraph)
    bl_idname = "LWSTUB"
    bl_label = "LW Stub Renderer"
    bl_use_alembic_procedural = False
    bl_use_custom_freestyle = False
    bl_use_eevee_viewport = True
    bl_use_gpu_context = False
    bl_use_image_save = False
    bl_use_materialx = False
    bl_use_postprocess = False
    bl_use_preview = False
    bl_use_shading_nodes_custom = False
    bl_use_spherical_stereo = False
    bl_use_stereo_viewport = False

    # Data to pass to and get from the stub renderer. Is probably very hacky global variable badness
    exported_objects: dict[bpy.types.ID, tuple[list[XMLNode], set[str]]] = None
    exported_object: XMLNode = None
    export_context: ExportContext = None
    export_success: bool = False

    def __init__(self):
        LoggingBlock.debug(f"[StubRenderer] Hello World!")

    def __del__(self):
        pass

    def render(self, depsgraph):
        # Preview is not supported
        if self.is_preview:
            LoggingBlock.warning(f"[StubRenderer] Tried to render preview!")
            return

        # We need an export context to do anything
        if LightwaveStubRenderEngine.export_context is None:
            LoggingBlock.warning(f"[StubRenderer] ExportContext was None!")
            return

        # Store render engine in context for setting frames
        LightwaveStubRenderEngine.export_context.engine = self
        export_success = False

        # Decide what to do
        if LightwaveStubRenderEngine.exported_objects is None:
            # No export map is given, so export camera
            LightwaveStubRenderEngine.exported_object = LightwaveStubRenderEngine.find_camera_instance_and_export(depsgraph, LightwaveStubRenderEngine.export_context.scene.camera)
        else:
            # If map is given, export objects
            LightwaveStubRenderEngine.export_layer_objects(depsgraph, LightwaveStubRenderEngine.export_context, LightwaveStubRenderEngine.exported_objects, depsgraph.view_layer)

        LightwaveStubRenderEngine.export_context.engine = None
        export_success = True
    
    @staticmethod
    def find_camera_instance_and_export(depsgraph: bpy.types.Depsgraph, camera: bpy.types.Object) -> XMLNode:
        for instance in depsgraph.object_instances:
            if (instance.object is not None) and instance.object.type == 'CAMERA' and instance.object.original == camera:
                return stub_export_camera(LightwaveStubRenderEngine.export_context, instance)
        
        LoggingBlock.error("No matching camera instance found")
        return None

    @staticmethod
    def export_layer_objects(depsgraph: bpy.types.Depsgraph, export_ctx: ExportContext, exported_objects: dict[bpy.types.ID, tuple[list[XMLNode], set[str]]], view_layer: bpy.types.ViewLayer) -> None:
        """Exports a single view layer by putting the not-yet-exported objects into the exported_objects dict"""

        if depsgraph.mode != 'RENDER':
            LoggingBlock.warning(f"Depsgraph was in mode '{depsgraph.mode}', but expected 'RENDER'")

        cur_frame = bpy.context.scene.frame_current
        start_frame = bpy.context.scene.frame_start

        # Export all object instances
        for instance in depsgraph.object_instances:
            object_eval = instance.object   # Get the underlying evaluated object for the instance

            if object_eval is None:
                continue
            if export_ctx.settings.use_selection and not object_eval.original.select_get():
                continue
            if not export_ctx.settings.use_selection and not instance.show_self:
                continue

            with LoggingBlock(f"exporting {'instance' if instance.is_instance else 'object'} with object name '{object_eval.name}'") as lb:
                if (object_eval.original in exported_objects) and (view_layer not in exported_objects[object_eval.original][1]):
                    # Object has already been exported, so just add this view layer to the set of view layers the object is present on
                    exported_objects[object_eval.original][1].add(view_layer)

                    lb.debug(f"Object is already exported for other view layers!")
                else:
                    # Object hasn't been exported yet, so we need to do that first before adding to map
                    obj_type = object_eval.type

                    res = None
                    
                    if obj_type in {'MESH', 'CURVE', 'SURFACE', 'META', 'FONT', 'CURVES'}:
                        frame_suffix = ""
                        
                        is_animated = object_eval.data.get("Animated", False)
                        if is_animated:
                            frame_suffix = "_frame_" + str(cur_frame)
                        
                        shapes: list[XMLNode] = export_ctx.registry.export(object_eval.original.data, lambda: export_shape(depsgraph, export_ctx, object_eval, frame_suffix))
                        meshcount: int = len(shapes)

                        if len(shapes) == 0:
                            export_ctx.warn(f"Entity {object_eval.name} has no material or shape and will be ignored")
                            
                        res = []
                            
                        for (mat_id, shape) in enumerate(shapes):
                            res.append(export_entity(export_ctx, instance, shape, mat_id))
                        
                    elif obj_type == "LIGHT" and export_ctx.settings.export_lights:
                        res = export_light(export_ctx, instance)

                    if res is None:
                        continue
                    
                    # Export particle system instances (hair) attached to object
                    hair_res = []
                    for ps in object_eval.particle_systems:
                        if ps.settings.type == 'HAIR' and ps.settings.render_type == 'PATH':
                            is_animated = ps.settings.get("Animated", False)
                            only_xml = (cur_frame > start_frame and is_animated)
                            frame_suffix = ""
                            if is_animated:
                                frame_suffix = "_frame_" + str(cur_frame)
                            lb.log(lb.LEVEL, f"Exporting hair particle system with name '{ps.settings.name}', adding .bob: {not only_xml}")
                            fur: XMLNode = export_fur(ps, export_ctx, export_ctx.registry.make_unique_id(ps.name + frame_suffix), only_xml)
                            hair_res.append(export_entity(export_ctx, instance, fur, ps.settings.material - 1))
                    

                    def _create_or_merge_list(a, b):
                        """Utility function to merge lists or single objects into a non-nested list"""
                        if isinstance(a, list):
                            if isinstance(b, list):
                                return a + b
                            else:
                                return a + [b]
                        else:
                            if isinstance(b, list):
                                return [a] + b
                            else:
                                return [a, b]

                    res = _create_or_merge_list(res, hair_res)

                    # Update exported_objects map with new exported nodes and view layers
                    if res is not None:
                        if object_eval.original in exported_objects:
                            exported_objects[object_eval.original] = (_create_or_merge_list(exported_objects[object_eval.original][0], res), exported_objects[object_eval.original][1])
                        else:
                            exported_objects[object_eval.original] = (res, {view_layer})

    def view_update(self, context, depsgraph):
        pass

    def view_draw(self, context, depsgraph):
        pass



particleCounter = 0
nonParticleCounter = 0

def export_technique(export_ctx: ExportContext) -> list[XMLNode]:
    if getattr(export_ctx.scene, 'cycles', None) is not None and bpy.context.engine == 'CYCLES':
        cycles = export_ctx.scene.cycles
        max_depth = cycles.max_bounces + 1
        clamp = max(cycles.sample_clamp_direct, cycles.sample_clamp_indirect)
        spp = cycles.samples
        use_denoising = cycles.use_denoising
    else:
        max_depth = 10
        clamp = 0
        spp = 64
        use_denoising = False
    
    name_prefix = f"{export_ctx.settings.scene_name_prefix.strip()}_" if len(export_ctx.settings.scene_name_prefix.strip()) > 0 else ""

    pt = XMLNode("integrator", type="pathtracer", depth=max_depth, samplingmode=export_ctx.settings.sampling_mode)
    pt.add("ref", id="scene")
    pt.add("image", id=f"{name_prefix}noisy")
    pt.add("sampler", type=export_ctx.settings.sampler_type.lower(), count=spp)

    if use_denoising:
        alb = XMLNode("integrator", type="aov", variable="albedo")
        alb.add("ref", id="scene")
        alb.add("image", id=f"{name_prefix}albedo")
        alb.add("sampler", type=export_ctx.settings.sampler_type.lower(), count=spp)

        norm = XMLNode("integrator", type="aov", variable="normals")
        norm.add("ref", id="scene")
        norm.add("image", id=f"{name_prefix}normals")
        norm.add("sampler", type=export_ctx.settings.sampler_type.lower(), count=spp)

        denoising = XMLNode("postprocess", type="denoising")
        denoising.add("ref", name="input", id=f"{name_prefix}noisy")
        denoising.add("ref", name="albedo", id=f"{name_prefix}albedo")
        denoising.add("ref", name="normals", id=f"{name_prefix}normals")
        denoising.add("image", id=f"{name_prefix}denoised")

        return [norm,alb,pt,denoising]
    
    return [pt]


def export_entity(export_ctx: ExportContext, inst: bpy.types.DepsgraphObjectInstance, shape: XMLNode, mat_id: int) -> XMLNode:
    """Exports the instance node for an object, including shape, transformation and material"""

    with LoggingBlock(f"exporting entity for instance of '{inst.object.name}'") as lb:
        instance_node = get_indexed_instance(export_ctx, inst)
        
        instance_node.add_child(shape)

        inst_mat = None
        if export_ctx.settings.export_materials:
            if mat_id < len(inst.object.material_slots):
                inst_mat = inst.object.material_slots[mat_id].material
                if inst_mat is not None:
                    instance_node.add_children(export_ctx.registry.export(inst_mat, lambda: export_material(export_ctx, inst_mat)))
                else:
                    export_ctx.warn(f"Obsolete material slot {mat_id} with instance {inst.object.data.name}. Maybe missing a material?")
                    instance_node.add_children(export_default_bsdf())
            else:
                #export_ctx.warn(f"Entity {inst.object.name} has no material")
                instance_node.add_children(export_default_bsdf())
        else:
            instance_node.add_children(export_default_bsdf())
            
        def _get_inner_transform(inst: bpy.types.DepsgraphObjectInstance):
            """Utility function passed to create_motion_blur_transform to get transformation matrix for instance"""
            trans_node = XMLNode("transform")
            trans_node.add("matrix", value=str_flat_matrix(inst.matrix_world))
            return trans_node

        # Motion Blur
        instance_node.add_child(create_motion_blur_transform(export_ctx, inst, _get_inner_transform))
        
        return instance_node

def get_indexed_instance(export_ctx: ExportContext, inst: bpy.types.DepsgraphObjectInstance) -> XMLNode:
    global particleCounter, nonParticleCounter
    if inst.particle_system:
        particleCounter += 1
        return XMLNode("instance", index = particleCounter)
    nonParticleCounter -= 1
    return XMLNode("instance", index = nonParticleCounter)

def export_view_layer(view_layer: bpy.types.ViewLayer):
    return XMLNode("viewlayer", name=re.sub("[^a-zA-Z0-9_\\- ]", "_", view_layer.name))

def export_single_render_layer(export_ctx: ExportContext, view_layer: bpy.types.ViewLayer) -> list[XMLNode]:
    """Exports a single render layer without adding view layer nodes"""

    with LoggingBlock(f"exporting single render layer {view_layer.name}") as lb:
        # Map that stores the XML nodes of already exported objects and the view layers they're visible on
        exported_objects: dict[bpy.types.ID, tuple[list[XMLNode], set[str]]] = {}

        # Set render engine to stub renderer and backup original render engine
        prev_engine = export_ctx.scene.render.engine
        export_ctx.scene.render.engine = 'LWSTUB'

        LightwaveStubRenderEngine.exported_objects = exported_objects
        LightwaveStubRenderEngine.export_context = export_ctx

        if export_ctx.scene.render.engine != 'LWSTUB':
            lb.error(f"Renderer is '{export_ctx.scene.render.engine}', expected 'LWSTUB'")
            return []
        
        bpy.ops.render.render(write_still=False, layer=view_layer.name, scene=export_ctx.scene.name, use_viewport=False)

        # Reset render engine
        export_ctx.scene.render.engine = prev_engine

        result = []

        for (conv, view_layers) in exported_objects.values():
            if isinstance(conv, list):
                result += [xml_node for xml_node in conv]
            else:
                result.append(conv)
        
        return result  

def export_render_layers(export_ctx: ExportContext) -> list[XMLNode]:
    """Exports the objects of add render layers using nodes for the view layers itself"""

    with LoggingBlock("exporting Render Layers") as lb:
        # Map that stores the XML nodes of already exported objects and the view layers they're visible on
        exported_objects: dict[bpy.types.ID, tuple[list[XMLNode], set[str]]] = {}

        # Set render engine to stub renderer and backup original render engine
        prev_engine = export_ctx.scene.render.engine
        export_ctx.scene.render.engine = 'LWSTUB'

        LightwaveStubRenderEngine.exported_objects = exported_objects
        LightwaveStubRenderEngine.export_context = export_ctx

        if export_ctx.scene.render.engine != 'LWSTUB':
            lb.error(f"Renderer is '{export_ctx.scene.render.engine}', expected 'LWSTUB'")
            return []

        result = []

        for view_layer in export_ctx.scene.view_layers:
            if not view_layer.use:
                lb.debug(f"Skipping view layer '{view_layer.name}' as it is disabled.")
                continue

            with LoggingBlock(f"exporting view layer '{view_layer.name}'") as lb:
                bpy.ops.render.render(write_still=False, layer=view_layer.name, scene=export_ctx.scene.name, use_viewport=False)
                result.append(export_ctx.registry.export(view_layer, lambda:export_view_layer(view_layer)))
        
        # Reset render engine
        export_ctx.scene.render.engine = prev_engine

        def _add_view_layers(xml_node: XMLNode, view_layers: set[bpy.types.ViewLayer]):
            """Utility function to add all view layers to object is visible on to it's instance node"""
            for view_layer in view_layers:
                xml_node.add_child(export_ctx.registry.export(view_layer, lambda:export_view_layer(view_layer)))
            return xml_node

        # Add view layer nodes to add  instances
        for (conv, view_layers) in exported_objects.values():
            if isinstance(conv, list):
                result += [_add_view_layers(xml_node, view_layers) for xml_node in conv]
            else:
                result += [_add_view_layers(conv, view_layers)]
        
        return result

def export_camera(export_ctx: ExportContext) -> XMLNode:
    """Exports the camera used to render the scene"""

    camera = export_ctx.scene.camera
    if camera is None:
        export_ctx.error("Your scene needs a camera!")
        return []
    
    # Set render engine to stub renderer and backup original render engine
    prev_engine = export_ctx.scene.render.engine
    export_ctx.scene.render.engine = 'LWSTUB'

    LightwaveStubRenderEngine.exported_objects = None
    LightwaveStubRenderEngine.export_context = export_ctx

    if export_ctx.scene.render.engine != 'LWSTUB':
        LoggingBlock.error(f"Renderer is '{export_ctx.scene.render.engine}', expected 'LWSTUB'")

    # Export using stub renderer to support motion blur
    bpy.ops.render.render(write_still=False, scene=export_ctx.scene.name, use_viewport=False)
    
    # Reset render engine
    export_ctx.scene.render.engine = prev_engine

    if LightwaveStubRenderEngine.exported_object is None:
        export_ctx.error("Couldn't export camera! Object was None")

    return [LightwaveStubRenderEngine.exported_object]

def export_scene(op, filepath, context, settings):
    with LoggingBlock("exporting Scene") as lb:
        # Root
        root = XMLRootNode()
        scene = XMLNode("scene", id="scene")

        global particleCounter, nonParticleCounter
        particleCounter = 0
        nonParticleCounter = 0

        # Create a path for meshes & textures
        rootPath = os.path.dirname(filepath)
        meshDir = os.path.join(rootPath, get_prefs().mesh_dir_name)
        texDir = os.path.join(rootPath, get_prefs().tex_dir_name)
        os.makedirs(meshDir, exist_ok=True)
        os.makedirs(texDir, exist_ok=True)

        registry = SceneRegistry()
        export_ctx = ExportContext(rootPath, settings, op, context.scene, registry)

        if export_ctx.scene.render.use_motion_blur and ('cycles' in export_ctx.scene) and (export_ctx.scene.cycles.motion_blur_position != 'START'):
            # We only support 'START' currently
            export_ctx.warn("Motion Blur position is not 'START', so export could look different")

        scene.add_children(export_camera(export_ctx))

        if settings.view_layer == 'ALL':
            scene.add_children(export_render_layers(export_ctx))
        else:
            vl_name = settings.view_layer.removeprefix("VL_")
            scene.add_children(export_single_render_layer(export_ctx, context.scene.view_layers[vl_name]))

        if settings.enable_background:
            scene.add_children(export_world_background(export_ctx))

            bpy.data.scenes["Scene"].render.use_compositing

        root.add_child(scene)
        
        if export_ctx.scene.render.use_compositing and export_ctx.scene.use_nodes:
            root.add_children(export_composition(export_ctx))
        else:
            # Leave simple export, if composition is not wanted/required
            root.add_children(export_technique(export_ctx))

        # Remove mesh & texture directory if empty
        try:
            if len(os.listdir(meshDir)) == 0:
                os.rmdir(meshDir)
            if len(os.listdir(texDir)) == 0:
                os.rmdir(texDir)
        except:
            pass  # Ignore any errors

        return root


def export_scene_to_file(op, filepath, context, settings):
    root = export_scene(op, filepath, context, settings)

    # Write the result into a file
    with open(filepath, 'w') as fp:
        fp.write(root.dump())
    
    logger.info(f"Finished exporting scene to file '{filepath}'")

    
def register():
    bpy.utils.register_class(LightwaveStubRenderEngine)

def unregister():
    bpy.utils.unregister_class(LightwaveStubRenderEngine)