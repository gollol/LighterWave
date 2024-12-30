# distutils: language = c++

from bob cimport *

cdef class PyVertex:
    cdef Vertex c_vert

    def __init__(self, p, uv, n):
        pos = Point(p[0], p[1], p[2])
        u = Vector2(uv[0], uv[1])
        no = Vector(n[0], n[1], n[2])
        self.c_vert = Vertex(pos, u, no)

cdef class PyMesh:
    cdef TriangleMesh c_mesh

    def addFace(self, i):
        indices = Vector3i(i[0], i[1], i[2])
        self.c_mesh.addFace(indices)
    
    def addVertex(self, PyVertex vs):
        self.c_mesh.addVertex(vs.c_vert)
    
    def save(self, path):
        self.c_mesh.saveMesh(path.encode('utf-8'))

cdef class PyHair:
    cdef Hair c_hair

    def addPoint(self, p):
        vec = Vector(p[0], p[1], p[2])
        self.c_hair.addPoint(vec)

cdef class PyFur:
    cdef Fur c_fur

    def __init__(self, numCurvePoints):
        self.c_fur = Fur(numCurvePoints)

    def addHair(self, PyHair h):
        self.c_fur.addHair(h.c_hair)
    
    def save(self, path):
        self.c_fur.saveFur(path.encode('utf-8'))