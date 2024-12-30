import bpy

from bpy.props import (
    StringProperty,
    BoolProperty
)
from bpy.types import AddonPreferences

package_name = __import__(__name__.split('.')[0])


class LightwavePreferences(AddonPreferences):
    bl_idname = package_name.__package__

    mesh_dir_name: StringProperty(
        name="Mesh Dir",
        description="Name for directory containing meshes",
        default="meshes",
    )
    tex_dir_name: StringProperty(
        name="Texture Dir",
        description="Name for directory containing textures",
        default="textures",
    )


    enable_debug_messages: BoolProperty(
        name="Enable Debug Messages",
        description="If true, print debug messages",
        default=True,
    )
    enable_log_file: BoolProperty(
        name="Enable Log File",
        description="If true, writes a log file with the name specified below",
        default=True,
    )
    log_file_name: StringProperty(
        name="Log File Name",
        description="File Name of the log file to generate",
        default="exporter.log",
    )

    def draw(self, context):
        layout = self.layout

        col = layout.column(heading="Directory Names")
        col.prop(self, "mesh_dir_name", text="Mesh")
        col.prop(self, "tex_dir_name", text="Texture")
        
        col = layout.column(heading="Logging")
        col.prop(self, "enable_debug_messages", text="Enable Debug Messages")
        col.prop(self, "enable_log_file", text="Enable Log File")
        col.prop(self, "log_file_name", text="Log File Name")


def register():
    bpy.utils.register_class(LightwavePreferences)


def unregister():
    bpy.utils.unregister_class(LightwavePreferences)


def get_prefs() -> LightwavePreferences:
    return bpy.context.preferences.addons[package_name.__package__].preferences
