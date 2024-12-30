import mathutils
import re
import logging
from datetime import datetime
import sys
import bpy

logger = logging.getLogger("lightwave_blender")

class LoggingBlock:
    """Context manager that provides the ability to indent log messages recursively based on task
    and profile the time taken for said task. Use this for easy logging."""
    indentation_level: int = 0
    logger: logging.Logger = logger
    
    SPACING_STR: str = "\t"     # Indentation string
    LEVEL: int = logging.DEBUG  # Log Level to use for profile messages
    ENABLED: bool = True        # Whether the LoggingBlock functionality is enabled in general
    
    def __init__(self, task_msg: str, active: bool = True):
        self.task_msg = task_msg
        self.active = active

        self.start_time = None
        self.finish_time = None
    
    def __enter__(self):
        if LoggingBlock.ENABLED and self.active:
            LoggingBlock.log(LoggingBlock.LEVEL, f"Started {self.task_msg}.")
            LoggingBlock.indentation_level += 1
            self.start_time = datetime.now()

        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        if LoggingBlock.ENABLED and self.active:
            self.finish_time = datetime.now()
            diff_secs = (self.finish_time - self.start_time).total_seconds()
            LoggingBlock.indentation_level -= 1

            if exc_type is not None:
                LoggingBlock.error(f"Error occured in logging block!", exc_info=(exc_type, exc_value, traceback))

            LoggingBlock.log(LoggingBlock.LEVEL, f"Finished {self.task_msg}! Took {round(diff_secs, 4)} seconds.")
    
    @classmethod
    def info(cls, msg, *args, **kwargs):
        cls.logger.info(msg, *args, **kwargs, extra={"indentation": cls.SPACING_STR * cls.indentation_level})
    
    @classmethod
    def warning(cls, msg, *args, **kwargs):
        cls.logger.warning(msg, *args, **kwargs, extra={"indentation": cls.SPACING_STR * cls.indentation_level})
    
    @classmethod
    def error(cls, msg, *args, **kwargs):
        cls.logger.error(msg, *args, **kwargs, extra={"indentation": cls.SPACING_STR * cls.indentation_level})
    
    @classmethod
    def critical(cls, msg, *args, **kwargs):
        cls.logger.critical(msg, *args, **kwargs, extra={"indentation": cls.SPACING_STR * cls.indentation_level})
    
    @classmethod
    def debug(cls, msg, *args, **kwargs):
        cls.logger.debug(msg, *args, **kwargs, extra={"indentation": cls.SPACING_STR * cls.indentation_level})
    
    @classmethod
    def log(cls, level, msg, *args, **kwargs):
        cls.logger.log(level, msg, *args, **kwargs, extra={"indentation": cls.SPACING_STR * cls.indentation_level})

def find_unique_name(used: set[str], name: str) -> str:
    unique_name = name
    index = 0

    while unique_name in used:
        unique_name = f"{name}.{index:03d}"
        index += 1
    
    used.add(unique_name)
    return unique_name

def get_image_name(view_layer: bpy.types.ViewLayer, pass_name: str, prefix: str = ""):
    vl_id = re.sub("[^a-zA-Z0-9_\\- ]", "_", view_layer.name)
    return f"{prefix}{vl_id}_{pass_name}"

def escape_identifier(name):
    return re.sub('[^a-zA-Z0-9_]', '_', name)


def flat_matrix(matrix):
    return [x for row in matrix for x in row]


def str_float(f: float):
    return "%.5g" % f


def str_flat_matrix(matrix):
    return ",  ".join([
        ",".join([ str_float(v) for v in row ])
        for row in matrix
    ])


def str_flat_array(array):
    if isinstance(array, float):
        return str_float(array)
    return ",".join([str_float(x) for x in array])


def orient_camera(matrix, skip_scale=False):
    # Y Up, Front -Z
    loc, rot, sca = matrix.decompose()
    return mathutils.Matrix.LocRotScale(loc, rot @ mathutils.Quaternion((0, 0, 1, 0)) @ mathutils.Quaternion((0, 0, 0, 1)), mathutils.Vector.Fill(3, 1) if skip_scale else sca)


def try_extract_node_value(value, default=0):
    try:
        return float(value)
    except:
        return default


def check_socket_if_constant(socket, value):
    if socket.is_linked:
        return False

    if socket.type == "RGBA" or socket.type == "VECTOR":
        return socket.default_value[0] == value and socket.default_value[1] == value and socket.default_value[2] == value
    else:
        return socket.default_value == value


def check_socket_if_black(socket):
    return check_socket_if_constant(socket, value=0)


def check_socket_if_white(socket):
    return check_socket_if_constant(socket, value=1)
