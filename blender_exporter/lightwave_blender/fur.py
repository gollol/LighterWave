import os

from .xml_node import XMLNode
from .utils import str_flat_array
from .addon_preferences import get_prefs
from .export_context import ExportContext
import importlib

def curves_save(filepath, ps):
    module = importlib.import_module('.bob', package='lightwave_blender')
    py_fur = getattr(module, 'PyFur', None)
    py_hair = getattr(module, 'PyHair', None)

    fur = py_fur(ps.settings.hair_step + 1)

    for h in ps.particles:
        hair = py_hair()
        for segment in h.hair_keys:
            hair.addPoint(segment.co)
        fur.addHair(hair)

    fur.save(filepath)

def export_fur(ps, export_ctx: ExportContext, name: str, only_xml: bool) -> XMLNode:
    scale = ps.settings.radius_scale
    root_radius = str(round(ps.settings.root_radius/2 * scale, 5))
    tip_radius = str(round(ps.settings.tip_radius/2 * scale, 5))

    rel_path = os.path.join(get_prefs().mesh_dir_name, f"{name}.bob")
    abs_path = os.path.join(export_ctx.path, rel_path)

    fur = XMLNode("shape", type="fur", rootRadius=root_radius, tipRadius=tip_radius, filename=rel_path.replace('\\', '/'))

    if not only_xml:    # only_xml -> file already exists (for animations)
        curves_save(abs_path, ps)
    
    return fur