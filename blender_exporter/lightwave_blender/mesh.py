# Export model from Blender and save to Binary Object Buffer (or .bob)

import bpy, importlib

def mesh_save(filepath, bm):
    module = importlib.import_module('.bob', package='lightwave_blender')
    mesh = getattr(module, 'PyMesh', None)
    vertex = getattr(module, 'PyVertex', None)

    uv_lay = bm.loops.layers.uv.active
    has_uv = uv_lay is not None

    mesh = mesh(has_uv)

    ply_vert_map = {}
    ply_vert_id = 0

    for f in bm.faces:
        pf = []

        for loop in f.loops:
            v = loop.vert

# Begin change diff to upstream Blender (6th. July 2022)
            uv = loop[uv_lay].uv[:] if has_uv else (0.0, 0.0)
            normal = v.normal if f.smooth else f.normal
            map_id = (v.co[:], normal[:], uv)

            # Identify unique vertex.
            if (_id := ply_vert_map.get(map_id)) is not None:
                pf.append(_id)
                continue
# End of change diff to upstream Blender

            mesh.addVertex(vertex(v.co, uv, normal))

            ply_vert_map[map_id] = ply_vert_id
            pf.append(ply_vert_id)

            ply_vert_id += 1

        mesh.addFace(pf)

    mesh.save(filepath)
