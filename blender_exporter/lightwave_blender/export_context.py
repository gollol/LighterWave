import bpy
from .registry import SceneRegistry
from .utils import LoggingBlock

class ExportContext(object):
    """Class to store all the data relevant for exporting the scene"""

    def __init__(self, path: str, settings, operator: bpy.types.Operator, scene: bpy.types.Scene, registry: SceneRegistry):
        self._path       = path      # The path that the scene is being exported to
        self._settings   = settings  # The object storing the export settings
        self._operator   = operator  # The export operator used for the export
        self._scene      = scene     # The scene to export
        self._registry   = registry  # The scene registry that stores already exported projects
        self.engine: bpy.types.RenderEngine = None # The render engine currently in use

    def debug(self, text: str):
        self._operator.report({'DEBUG'}, text)
        LoggingBlock.debug(text)

    def warn(self, text: str):
        self._operator.report({'WARNING'}, text)
        LoggingBlock.warning(text)

    def error(self, text: str):
        self._operator.report({'ERROR'}, text)
        LoggingBlock.error(text)

    def info(self, text: str):
        self._operator.report({'INFO'}, text)
        LoggingBlock.info(text)

    #
    #   Properties to be accessed from the outside
    #

    @property
    def path(self) -> str:
        """The path that the scene is being exported to"""
        return self._path

    @property
    def settings(self):
        """The object storing the export settings"""
        return self._settings

    @property
    def operator(self) -> bpy.types.Operator:
        """The export operator used for the export"""
        return self._operator

    @property
    def scene(self) -> bpy.types.Scene:
        """The scene to export"""
        return self._scene

    @property
    def registry(self) -> SceneRegistry:
        """The scene registry that stores already exported projects"""
        return self._registry
    
    #
    #   Utility functions
    #

    def is_view_layer_selected(self, view_layer: bpy.types.ViewLayer):
        if self._settings.view_layer == 'ALL':
            return True

        vl_name = self._settings.view_layer.removeprefix("VL_")

        return self._scene.view_layers[vl_name] == view_layer