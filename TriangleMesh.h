#ifndef POLYSLICER_TRIANGLEMESH_H
#define POLYSLICER_TRIANGLEMESH_H

class TriangleMesh {

public:

    TriangleMesh() : bottomLeftVertex(999999,999999,999999), topRightVertex(-999999,-999999,-999999) { meshSize = 0;}

    size_t size() const {
        return meshSize;
    }

    void push_back(Triangle &t) {
        meshSize++;
        triangles.push_back(t);
        for (size_t i = 0; i < 3; ++i) {
            if (t.v[i].x < bottomLeftVertex.x) { bottomLeftVertex.x = t.v[i].x; }
            if (t.v[i].y < bottomLeftVertex.y) { bottomLeftVertex.y = t.v[i].y; }
            if (t.v[i].z < bottomLeftVertex.z) { bottomLeftVertex.z = t.v[i].z; }
            if (t.v[i].x > topRightVertex.x) { topRightVertex.x = t.v[i].x; }
            if (t.v[i].y > topRightVertex.y) { topRightVertex.y = t.v[i].y; }
            if (t.v[i].z > topRightVertex.z) { topRightVertex.z = t.v[i].z; }
        }
    }

    Vector3 meshAABBSize() const {
        return Vector3 ( topRightVertex.x - bottomLeftVertex.x,
                         topRightVertex.y - bottomLeftVertex.y,
                         topRightVertex.z - bottomLeftVertex.z );
    }

    const std::vector<Triangle>& getTriangles() const { return triangles; }

    Vector3 getBottomLeftVertex() const { return bottomLeftVertex; }

    Vector3 getTopRightVertex() const { return topRightVertex; }

public:

    int meshSize;
    Vector3 bottomLeftVertex;
    Vector3 topRightVertex;
    std::vector<Triangle> triangles;
};

#endif //POLYSLICER_TRIANGLEMESH_H
