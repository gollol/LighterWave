import bpy

from bpy.props import (
    BoolProperty,
    StringProperty,
    EnumProperty
)
from bpy_extras.io_utils import (
    ExportHelper
)
from bpy_extras.wm_utils.progress_report import (
    ProgressReport
)
from collections import namedtuple

from .addon_preferences import get_prefs

from .exporter import *

import logging, os

from .utils import logger


class ExportLightwave(bpy.types.Operator, ExportHelper):
    """Export scene to Lightwave"""

    bl_idname = "export_scene.lightwave"
    bl_label = "Export Lightwave Scene"
    bl_description = "Export scene to Lightwave"
    bl_options = {'PRESET'}

    filename_ext = ".xml"
    filter_glob: StringProperty(
        default="*.xml",
        options={'HIDDEN'}
    )

    use_selection: BoolProperty(
        name="Selection Only",
        description="Export selected objects only",
        default=False,
    )

    animations: BoolProperty(
        name="Export Animations",
        description="If true, writes .xml for each frame in the animation. Mark each animated object via a custom property 'Animated' set to true!",
        default=False,
    )

    export_materials: BoolProperty(
        name="Export Materials",
        description="If true, materials will be exported, else a single diffuse material will be used",
        default=True,
    )

    export_lights: BoolProperty(
        name="Export Lights",
        description="If true, lights will be exported",
        default=True,
    )

    enable_background: BoolProperty(
        name="Export Background",
        description="If true, background will be exported as a light",
        default=True,
    )

    enable_camera: BoolProperty(
        name="Export Camera",
        description="If true, active camera will be exported",
        default=True,
    )

    enable_integrator: BoolProperty(
        name="Export Integrator",
        description="If true, current integrator technique will be mapped to Lightwave",
        default=True,
    )

    copy_images: BoolProperty(
        name="Copy all Images",
        description="If true, copy all images next to the scene file, not only Generated or Packed images",
        default=True,
    )

    enable_area_lights: BoolProperty(
        name="Export light as area",
        description="If true, lights will be exported as area lights. If your renderer does not support area lights consider turning this off",
        default=True,
    )

    overwrite_existing_meshes: BoolProperty(
        name="Overwrite existing meshes",
        description="Disable this to speed up export if you have not modified any geometry",
        default=True,
    )

    overwrite_existing_textures: BoolProperty(
        name="Overwrite existing textures",
        description="Disable this to speed up export if you have not modified any textures",
        default=True,
    )

    sampling_mode: EnumProperty(
        name="Sampling Mode",
        description="Decides which sampling mode will be used in the exported lightwave scene file",
        items=[
            ("BSDF", "BSDF", "Use only the BSDFs of objects for rendering the image"),
            ("NEE", "NEE", "Additionally use Next Event Estimation next to default BSDF sampling"),
            ("MIS", "MIS", "Use Multiple Importance Sampling to combine default BSDF and NEE sampling"),
        ],
        default="MIS"
    )

    sampler_type: EnumProperty(
        name="Sampler Type",
        description="Decides which sampler to use in the exported lighwave scene file",
        items=[
            ("INDEPENDENT", "Independent", "Use the independent sampler"),
            ("HALTON", "Halton", "Use the halton sampler")
        ],
        default="INDEPENDENT"
    )

    environment_map_sampling_mode: EnumProperty(
        name="Environment Map Sampling Mode",
        description="Decides the technique that is used to sample the environment map!",
        items=[
            ("SIMPLE", "SIMPLE", "The environment map is sampled uniformly"),
            ("STANDARD", "STANDARD", "The environment map is sampled with the standard sampler"),
            ("WARP", "WARP", "The environment map is sampled with the warp sampling"),
        ],
        default="STANDARD"
    )

    scene_name_prefix: StringProperty(
        name="Name Prefix",
        description="If not empty, sets the scene name prefix for the names of the image files generated from the lighwave scene file",
        default=""
    )

    def get_view_layer_items(self, context) -> list:
        """Get list of view layers dynamically"""
        if context is None:
            return []
        
        items = [
            ("ALL", "All Layers", "Export all (activated) view layers", 0),
            None,
            (f"VL_{context.view_layer.name}", context.view_layer.name, f"Export only layer '{context.view_layer.name}'", 1)
        ]

        for view_layer in context.scene.view_layers:
            if not view_layer.use or (view_layer == context.view_layer):
                continue
            
            items.append(
                (f"VL_{view_layer.name}", view_layer.name, f"Export only layer '{view_layer.name}'")
            )
        
        return items

    view_layer: EnumProperty(
        name="View/Render Layer",
        description="""Which view/render layer(s) should be exported to the exported lighwave scene file.
In the case a single layer is selected, no view layer nodes will be added to instances.""",
        items=get_view_layer_items,
        default=1
    )

    check_extension = True

    def execute(self, context):
        keywords = self.as_keywords(
            ignore=(
                "filepath",
                "filter_glob",
                "animations",
                "check_existing"
            ),
        )

        # Setup logging
        logger.setLevel(logging.DEBUG)
        logger.propagate = False

        # We could still have handlers, so remove them to avoid duplicate messages
        for h in logger.handlers:
            logger.removeHandler(h)

        formatter = logging.Formatter(
            fmt="[%(asctime)s] [%(name)s/%(levelname)-5.5s]%(indentation)s %(message)s",
            datefmt="%H:%M:%S",
            defaults={"indentation": ""}
        )
        consoleHandler = logging.StreamHandler(sys.stdout)
        consoleHandler.setFormatter(formatter)
        if not get_prefs().enable_debug_messages:
            consoleHandler.setLevel(logging.INFO)
        logger.addHandler(consoleHandler)

        if get_prefs().enable_log_file:
            logfileHandler = logging.FileHandler(filename=os.path.join(os.path.dirname(self.filepath), get_prefs().log_file_name), mode='w')
            logfileHandler.setFormatter(formatter)
            logfileHandler.setLevel(logging.DEBUG)
            logger.addHandler(logfileHandler)

        # Exit edit mode before exporting, so current object states are exported properly.
        if bpy.ops.object.mode_set.poll():
            bpy.ops.object.mode_set(mode='OBJECT')

        settings = namedtuple("Settings", keywords.keys())(*keywords.values())

        try:
            with ProgressReport(context.window_manager) as progress:
                if self.animations is True:
                    scene_frames = range(
                        context.scene.frame_start, context.scene.frame_end + 1)
                    progress.enter_substeps(len(scene_frames))
                    for frame in scene_frames:
                        context.scene.frame_set(frame)
                        progress.enter_substeps(1)
                        export_scene_to_file(self, self.filepath.replace(
                            '.xml', f'{frame:04}.xml'), context, settings=settings)
                    progress.leave_substeps()
                else:
                    export_scene_to_file(self, self.filepath, context, settings=settings)
            
            self.report({'INFO'}, "Exported scene to %s" % (self.filepath))
        finally:
            # We're done, so remove handlers again. This also unlocks the exporter log file
            for h in logger.handlers:
                logger.removeHandler(h)

        return {'FINISHED'}

    def draw(self, context):
        pass


class LIGHTWAVE_PT_export_include(bpy.types.Panel):
    bl_space_type = 'FILE_BROWSER'
    bl_region_type = 'TOOL_PROPS'
    bl_label = "Include"
    bl_parent_id = "FILE_PT_operator"

    @classmethod
    def poll(cls, context):
        sfile = context.space_data
        operator = sfile.active_operator

        return operator.bl_idname == "EXPORT_SCENE_OT_lightwave"

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False  # No animation.

        sfile = context.space_data
        operator = sfile.active_operator

        col = layout.column(heading="Limit to")
        col.prop(operator, 'use_selection')

        layout.separator()
        col = layout.column(heading="Export")
        col.prop(operator, 'animations', text="Animations")
        col.prop(operator, 'export_materials', text="Materials")
        col.prop(operator, 'export_lights', text="Lights")
        col.prop(operator, 'enable_background', text="Background")
        col.prop(operator, 'enable_camera', text="Camera")
        col.prop(operator, 'enable_integrator', text="Integrator")
        col.prop(operator, 'view_layer', text="Render Layer")

        layout.separator()
        col = layout.column(heading="Images")
        col.prop(operator, 'copy_images')

        layout.separator()
        col = layout.column(heading="Overwrite")
        col.prop(operator, 'overwrite_existing_meshes', text="Existing Meshes")
        col.prop(operator, 'overwrite_existing_textures', text="Existing Textures")

        layout.separator()
        col = layout.column(heading="Features")
        col.prop(operator, 'enable_area_lights')
        col.prop(operator, 'sampling_mode')
        col.prop(operator, 'sampler_type')
        col.prop(operator, 'environment_map_sampling_mode')

        layout.separator()
        col = layout.column()
        col.prop(operator, 'scene_name_prefix')



def menu_func_export(self, context):
    self.layout.operator(ExportLightwave.bl_idname, text="Lightwave (.xml)")


classes = (
    ExportLightwave,
    LIGHTWAVE_PT_export_include
)


def register():
    for cls in classes:
        bpy.utils.register_class(cls)
    bpy.types.TOPBAR_MT_file_export.append(menu_func_export)


def unregister():
    bpy.types.TOPBAR_MT_file_export.remove(menu_func_export)
    for cls in classes:
        bpy.utils.unregister_class(cls)


if __name__ == "__main__":
    register()
