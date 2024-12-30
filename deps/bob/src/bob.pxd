cdef extern from "math.hpp" namespace "internal":
    cdef cppclass Point:
        Point() except +
        Point(float x, float y, float z) except +
    
    cdef cppclass Vector2:
        Vector2() except +
        Vector2(float x, float y) except +
    
    cdef cppclass Vector:
        Vector() except +
        Vector(float x, float y, float z) except +
    
    cdef cppclass Vector3i:
        Vector3i() except +
        Vector3i(int x, int y, int z) except +

cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        pass

cdef extern from "mesh.cpp" namespace "internal":
    cdef cppclass TriangleMesh:
        TriangleMesh() except +
        TriangleMesh(uvs) except +

        void addVertex(Vertex)
        void addFace(Vector3i)
    
        void saveMesh(string path)
    
    cdef cppclass Vertex:
        Vertex() except +
        Vertex(Point position, Vector2 uv, Vector normal) except +

cdef extern from "utils.hpp" namespace "internal":    
    cdef cppclass Hair:
        Hair() except +
        void addPoint(Vector point)

    cdef cppclass Fur:
        Fur() except +
        Fur(int keyPoints) except +

        void addHair(Hair hair)
        void saveFur(string path)