import math
from .utils import *
from .xml_node import XMLNode
from .export_context import ExportContext
import bpy
from .transform import create_motion_blur_transform

def stub_export_camera(export_ctx: ExportContext, instance: bpy.types.DepsgraphObjectInstance) -> XMLNode:
    with LoggingBlock("exporting camera") as lb:
        camera = export_ctx.scene.camera
        render = export_ctx.scene.render
        res_x = int(render.resolution_x * render.resolution_percentage * 0.01)
        res_y = int(render.resolution_y * render.resolution_percentage * 0.01)

        # TODO: Other types?
        camera_node = XMLNode("camera", type="perspective")

        camera_node.add("integer", name="width", value=res_x)
        camera_node.add("integer", name="height", value=res_y)

        camera_node.add("float", name="fov", value=math.degrees(
            2 * math.atan(camera.data.sensor_width / (2 * camera.data.lens))))
        camera_node.add("string", name="fovAxis", value="x" if render.resolution_x > render.resolution_y else "y")

        # camera_node.add("float", name="nearClip", value=camera.data.clip_start)
        # camera_node.add("float", name="farClip", value=camera.data.clip_end)

        def _get_inner_transform(inst: bpy.types.DepsgraphObjectInstance):
            trans_node = XMLNode("transform")
            trans_node.add("matrix", value=str_flat_matrix(orient_camera(instance.matrix_world, skip_scale=True)))
            return trans_node

        # Motion Blur
        camera_node.add_child(create_motion_blur_transform(export_ctx, instance, _get_inner_transform))

        return camera_node
