#ifndef POLYSLICER_TRIANGLEMESH_H
#define POLYSLICER_TRIANGLEMESH_H

#include <cstddef>
#include <vector>

#include "Triangle.h"


class TriangleMesh {

public:
    float roundToNearestEpsilon(float x, double epsilon, int mod, int rem) {
        double y = round((double) x / (mod * epsilon));
        double z = (y * mod + rem) * epsilon;
        return (float) z;
    }

    Vector3 roundVector3(float x, float y, float z, double epsilon) {
        Vector3 p;
        p.x = roundToNearestEpsilon(x, epsilon, 2, 0);
        p.y = roundToNearestEpsilon(y, epsilon, 2, 0);
        p.z = roundToNearestEpsilon(z, epsilon, 2, 0);
        return p;
    }

    Triangle makeRoundedTriangle(
            float n0, float n1, float n2,
            float f0, float f1, float f2,
            float f3, float f4, float f5,
            float f6, float f7, float f8,
            double epsilon) {

        return Triangle(Vector3(n0, n1, n2),
                        roundVector3(f0, f1, f2, epsilon),
                        roundVector3(f3, f4, f5, epsilon),
                        roundVector3(f6, f7, f8, epsilon));
    }

    TriangleMesh() : bottomLeftVertex(999999,999999,999999), topRightVertex(-999999,-999999,-999999) { meshSize = 0;}

    TriangleMesh(const char *stlFile, double eps) {

        int numberDegenerateTriangles = 0;

        FILE *f = fopen(stlFile, "rb");
        if (!f) {
            throw std::invalid_argument("File: " + std::string(stlFile) + "  - Does not exist!");
        }
        char header[80];
        int numberTriangles;
        int err;
        err = fread(header, 80, 1, f);
        err = fread((void *) &numberTriangles, 4, 1, f);
        float v[12]; /* normal = 3, vertices = 9  */
        unsigned short uint16;

        for (size_t i = 0; i < numberTriangles; ++i) {
            for (size_t j = 0; j < 12; ++j) {
                err = fread((void *) &v[j], sizeof(float), 1, f);
            }
            err = fread((void *) &uint16, sizeof(unsigned short), 1, f); // spacer between successive faces
            Triangle triangle = makeRoundedTriangle(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11],
                                              eps);
            if (!triangle.isDegenerate()) {
                push_back(triangle);
            } else {
                numberDegenerateTriangles++;
            }
        }
        fclose(f);

        std::cout << "Number triangles: " << numberTriangles << std::endl;
        std::cout << "Number degenerate triangles: " << numberDegenerateTriangles << std::endl;
    }

    size_t size() const {
        return meshSize;
    }

    void push_back(Triangle &triangle) {
        meshSize++;
        triangles.push_back(triangle);
        for (size_t i = 0; i < 3; ++i) {
            if (triangle.v[i].x < bottomLeftVertex.x) { bottomLeftVertex.x = triangle.v[i].x; }
            if (triangle.v[i].y < bottomLeftVertex.y) { bottomLeftVertex.y = triangle.v[i].y; }
            if (triangle.v[i].z < bottomLeftVertex.z) { bottomLeftVertex.z = triangle.v[i].z; }
            if (triangle.v[i].x > topRightVertex.x) { topRightVertex.x = triangle.v[i].x; }
            if (triangle.v[i].y > topRightVertex.y) { topRightVertex.y = triangle.v[i].y; }
            if (triangle.v[i].z > topRightVertex.z) { topRightVertex.z = triangle.v[i].z; }
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
