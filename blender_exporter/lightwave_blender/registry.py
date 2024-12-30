import bpy
import os
import re
from collections.abc import Callable

from .utils import find_unique_name, logger, LoggingBlock
from .xml_node import XMLNode
class Token(object):
    """Class that is used to identify an object in the scene registry, which doesn't have an original property to compare"""

    def __init__(self, object, **data):
        self.object = object
        self.data = data
    
    def get_id_name(self):
        """gets a string id name of the object and the values of the associated data"""
        id_name: str = str(self.object)

        for _, value in self.data.items():
            id_name = id_name + "_" + str(value)
        
        return id_name

    def __eq__(self, other):
        if isinstance(other, self.__class__) and (self.data.keys() == other.data.keys()):
            res = (self.object == other.object)
            for key, value in self.data.items():
                res = res and (value == other.data[key])
            return res
        else:
            return False
    
    def __hash__(self):
        h = hash(self.object)
        for key, value in self.data.items():
            h = hash((h, key, value))
        return h
    
    def __str__(self):
        return self.get_id_name()
    
    @property
    def name(self):
        return self.get_id_name()


class SceneRegistry(object):
    """Class that stores already exported objects (original objects) and their corresponding exported XML Nodes"""

    def __init__(self):
        # The map of original objects to the exported XML Node
        self.converted: dict[bpy.types.ID | Token, XMLNode | list[XMLNode]] = {}
        # The set of used ID names
        self.used_ids: set[str] = set()
    
    def make_unique_id(self, id_name: str) -> str:
        """Creates a unique ID (relative to this registry) from a given ID name and the already used IDs stored in the registry"""

        id_name = re.sub("[^a-zA-Z0-9_\\- ]", "_", id_name)
        return find_unique_name(self.used_ids, id_name)

    def _build_reference(self, node: XMLNode, id_name: str) -> XMLNode:
        """Builds a reference to the given node and sets it's id attribute"""
        

        if "id" not in node.attributes:
            node.attributes["id"] = self.make_unique_id(id_name)

        kw = {}
        if "name" in node.attributes:
            kw["name"] = node.attributes["name"]
        kw["id"] = node.attributes["id"]

        return XMLNode("ref", **kw)

    def export(self, original: bpy.types.ID | Token, export_fn: Callable[[], XMLNode | list[XMLNode]]) -> XMLNode | list[XMLNode]:
        """Exports an object using the provided export function or,
        if already present in the converted map, returns the existing exported node"""

        if original in self.converted:
            LoggingBlock.debug(f"Original '{original.name}' exists in the registry! Building reference...")
            # Object was already exported, so build a reference to existing one
            conv = self.converted[original]

            if isinstance(conv, list):
                return [ self._build_reference(node, original.name) for node in conv ]
            return self._build_reference(conv, original.name)

        # If not already exported, export the object and add to map
        self.converted[original] = export_fn()
        LoggingBlock.debug(f"Original '{original.name}' was not in registry yet and has been added")

        return self.converted[original]
    
    def dump_log(self):
        """Debugging function to dump all entries stored in registry to the log"""
        
        with LoggingBlock("dumping scene registry contents") as lb:
            lb.debug("Used IDs:")
            for uid in self.used_ids:
                lb.debug(f"  {uid}")
            lb.debug("")

            lb.debug("Converted Elements:")

            for obj, conv in self.converted.items():
                lb.debug(f"  {type(obj).__name__}: {str(obj)} => {(str([n.dump().replace(chr(10), ' ') for n in conv])) if isinstance(conv, list) else conv.dump().replace(chr(10), ' ')}")