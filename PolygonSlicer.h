#ifndef POLYSLICER_POLYGONSLICER_H
#define POLYSLICER_POLYGONSLICER_H

#include "TriangleMesh.h"
#include "Triangle.h"
#include "LineSegment.h"
#include "Contour.h"
#include "BoundingBox.h"
#include "Layer.h"


class PolygonSlicer {
public:
    std::vector<Layer *> sliceModel(TriangleMesh &mesh, double thickness, double epsilon = 0.0001);

private:
    typedef struct node {
        Triangle t;
        struct node *next;
        struct node *prev;
    } Mesh_Triangle_Node_t;

    typedef struct _list {
        Mesh_Triangle_Node_t *head;
        Mesh_Triangle_Node_t *tail;
    } Mesh_Triangle_List_t;

    Mesh_Triangle_List_t* Mesh_Triangle_List_create (void);

    void Mesh_Triangle_List_insert (Triangle t, Mesh_Triangle_List_t *L);

    void Mesh_Triangle_List_union (Mesh_Triangle_List_t *L1, Mesh_Triangle_List_t *L2);

    Mesh_Triangle_Node_t* Mesh_Triangle_List_remove (Mesh_Triangle_List_t *L, Mesh_Triangle_Node_t *node);

    std::vector<Layer *> TrivialSlicing(const TriangleMesh &mesh, std::vector<float> &planes);

    std::vector<Layer *> IncrementalSlicing(const TriangleMesh& mesh, std::vector<float> &planes, double delta);

    void TrivialLoopClosure(std::vector<LineSegment> lineSegments, std::vector<Contour> &contour);

    Vector3 TrivialClosure_find_A(std::vector<LineSegment> &segs, Vector3 u);

    LineSegment R3_Mesh_Triangle_slice(Triangle t, float Z);

    Vector3 R3_Mesh_Side_slice(Vector3 vi, Vector3 vj, float Z);

    LineSegment R3_Mesh_Triangle_slice (Mesh_Triangle_Node_t *t, float Z);

    std::vector<float> computePlanes(TriangleMesh &mesh, double thickness, double epsilon);

    Mesh_Triangle_List_t** IncrementalSlicing_buildLists (const TriangleMesh& mesh, std::vector<float> P, double delta);

    Triangle makeRoundedTriangle(
            float n0, float n1, float n2,
            float f0, float f1, float f2,
            float f3, float f4, float f5,
            float f6, float f7, float f8,
            double epsilon);

    Vector3 roundVector3(float x, float y, float z, double epsilon);

    float roundToNearestEpsilon(float x, double epsilon, int mod, int rem);

    LineSegment create_ray(Vector3 point, BoundingBox bbox, int index);

    bool is_inside(LineSegment line, Vector3 point);

    bool ray_intersect(LineSegment ray, LineSegment side);

    void update_bounding_box(Vector3 point, BoundingBox *bbox);

    bool insided_bounding_box(Vector3 point, BoundingBox bb);

    bool contains(Vector3 point, BoundingBox bb, std::vector<LineSegment> sides, int index);

    void add_point(Vector3 p1, Vector3 p2, std::vector<LineSegment> &t, BoundingBox *bbox, bool first, int index);

    void ray_casting(std::vector<Contour> &polygons);

    };


#endif //POLYSLICER_POLYGONSLICER_H
