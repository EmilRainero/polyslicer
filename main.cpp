#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include "Timer.h"
#include "Vector3.h"
#include "Triangle.h"
#include "TriangleMesh.h"
#include "LineSegment.h"


double eps = 0.0001; // round all vertices in stl file to within this epsilon



float roundToNearestEpsilon(float x, double epsilon, int mod, int rem) {
    double y = round((double) x / (mod * epsilon));
    double z = (y * mod + rem) * epsilon;
    return (float) z;
}

Vector3 roundVector3(float x, float y, float z, double eps) {
    Vector3 p;
    p.x = roundToNearestEpsilon(x, eps, 2, 0);
    p.y = roundToNearestEpsilon(y, eps, 2, 0);
    p.z = roundToNearestEpsilon(z, eps, 2, 0);
    return p;
}

Triangle makeRoundedTriangled(
        float n0, float n1, float n2,
        float f0, float f1, float f2,
        float f3, float f4, float f5,
        float f6, float f7, float f8,
        double eps) {

    return Triangle(Vector3(n0, n1, n2),
                    roundVector3(f0, f1, f2, eps),
                    roundVector3(f3, f4, f5, eps),
                    roundVector3(f6, f7, f8, eps));
}

int stlToMeshInMemory(const char *stlFile, TriangleMesh &mesh) {

    int numberDegenerateTriangles = 0;

    FILE *f = fopen(stlFile, "rb");
    if (!f) {
        return -1;
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
        Triangle t = makeRoundedTriangled(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], eps);
        if (!t.degenerate()) {
            mesh.push_back(t);
        } else {
            numberDegenerateTriangles++;
        }
    }
    fclose(f);

    std::cout << "Number triangles: " << numberTriangles << std::endl;
    std::cout << "Number degenerate triangles: " << numberDegenerateTriangles << std::endl;
    return 0;
}

std::vector<float> compute_planes(TriangleMesh &mesh, double thickness) {
    std::vector<float> planes;

    double model_zmax = mesh.getTopRightVertex().z;
    double model_zmin = mesh.getBottomLeftVertex().z;

    double rounded_thickness = roundToNearestEpsilon(thickness, eps, 2, 0); /*Plane spacing even multiple of {eps}*/
    double P0 = roundToNearestEpsilon(model_zmin - rounded_thickness, eps, 2, 1); /*First plane odd multiple of {eps}.*/

    int no_planes = 1 + (int) ((model_zmax + rounded_thickness - P0) / rounded_thickness);

    std::cout << "Model Z: [" << model_zmin << "  " << model_zmax << "] = " << model_zmax - model_zmin << " height"
              << std::endl;

    for (size_t i = 0; i < no_planes; i++) {
        float Pi = (float) (P0 + i * rounded_thickness);
        if ((Pi > model_zmin) && (Pi < model_zmax)) {
            planes.push_back(Pi);
        }
    }

    return planes;
}

typedef struct node {
    Triangle t;
    struct node *next;
    struct node *prev;
} Mesh_Triangle_Node_t;

typedef struct _list {
    Mesh_Triangle_Node_t *head;
    Mesh_Triangle_Node_t *tail;
} Mesh_Triangle_List_t;

Vector3 R3_Mesh_Side_slice(Vector3 vi, Vector3 vj, float Z) {
    double dx = vj.x - vi.x;
    double dy = vj.y - vi.y;
    double dz = vj.z - vi.z;
//    assert(dz != 0);
    double frac = (Z - vi.z) / dz;
    float xint = (float) (frac * dx + (double) vi.x);
    float yint = (float) (frac * dy + (double) vi.y);
    return (Vector3) {.x = xint, .y = yint, .z = Z};
}

LineSegment R3_Mesh_Triangle_slice(Triangle t, float Z) {
//    assert((t.zMin < Z) && (t.zMax > Z));
    int np = 0; /* Number of segment endpoints found */
    LineSegment segment;
    for (int i = 0; i < 3; i++) {
        /* Get side {i} of triangle: */
        int j = (i == 2 ? 0 : i + 1);
        Vector3 vi = (t.v[i]);
        Vector3 vj = (t.v[j]);
        /* Check for intersection of plane with {vi--vj}. */
        /* Must consider segment closed at bottom and open at top in case {Z} goes through a vertex. */
        float vzMin = (vi.z < vj.z ? vi.z : vj.z);
        float vzMax = (vi.z > vj.z ? vi.z : vj.z);
        if ((vzMin <= Z) && (vzMax > Z)) {
            Vector3 p = R3_Mesh_Side_slice(vi, vj, Z);
//            assert(np < 2);
            segment.v[np] = p;
            np++;
        }
    }
//    assert(np == 2);
    return segment;
}

Vector3 TrivialClosure_find_A(std::vector<LineSegment> &segs, Vector3 u) {
    for (std::vector<LineSegment>::iterator j = segs.begin(); j != segs.end(); ++j) {
        Vector3 p0 = (*j).v[0];
        Vector3 p1 = (*j).v[1];
        if ((p0.x == u.x) && (p0.y == u.y)) {
            segs.erase(j);
            return p1;
        } else if ((p1.x == u.x) && (p1.y == u.y)) {
            segs.erase(j);
            return p0;
        }
    }
    return {+INFINITY, +INFINITY, +INFINITY};
}

typedef struct _bounding_box {
    double xMin;
    double xMax;
    double yMin;
    double yMax;
    bool firstPoint;
} bounding_box;

typedef struct _contour {
    bool external;
    bool clockwise;
    std::vector<Vector3> P;
} contour;

void TrivialLoopClosure(std::vector<LineSegment> segs, std::vector<contour> &contours) {

    while (segs.size() > 0) {

        /* Get another contour: */
        std::vector<Vector3> P;

        /* Get the first segment: */
        Vector3 last;
        Vector3 current;
        Vector3 prev;
        {
            std::vector<LineSegment>::iterator i = segs.begin();
            last = (*i).v[0];
            current = (*i).v[1];
            {
                P.push_back(last);
            }
            prev = last;
            segs.erase(i);
        }

        /* Get additional segments until loop is closed: */
        bool open = false;
        do {
            /* Find and delete another segment with endpoint {current}, advance {prev,current} */
            Vector3 next = TrivialClosure_find_A(segs, current);
            //v3 next = TrivialClosure_find_B (segs, current, plane);
            if ((next.x == +INFINITY) || (next.y == +INFINITY)) {
                /* Open contour?! */
                open = true;
                break;
            }
            {
                P.push_back(current);
            }
            prev = current;
            current = next;
        } while ((current.x != last.x) || (current.y != last.y));
        if (!open) {
            P.push_back(last);
        }
        contours.push_back({false, false, P});
    }
}

LineSegment create_ray(Vector3 point, bounding_box bb, int index) {
    /* Create outside point: */
    double epsilon = (bb.xMax - bb.xMin) / 100.0;
    Vector3 outside(bb.xMin - epsilon, bb.yMin);
    LineSegment v(outside, point, index);
    return v;
}

bool is_inside(LineSegment line, Vector3 point) {
    double maxX = (line.v[0].x > line.v[1].x) ? line.v[0].x : line.v[1].x;
    double minX = (line.v[0].x < line.v[1].x) ? line.v[0].x : line.v[1].x;
    double maxY = (line.v[0].y > line.v[1].y) ? line.v[0].y : line.v[1].y;
    double minY = (line.v[0].y < line.v[1].y) ? line.v[0].y : line.v[1].y;
    if ((point.x >= minX && point.x <= maxX) && (point.y >= minY && point.y <= maxY)) {
        return true;
    }
    return false;
}

bool ray_intersect(LineSegment ray, LineSegment side) {
    Vector3 intersectPoint;
    /* If both vectors aren't from the kind of x=1 lines then go into: */
    if (!ray.vertical && !side.vertical) {
        /* Check if both vectors are parallel. If they are parallel then no intersection point will exist: */
        if (ray.a - side.a == 0) {
            return false;
        }
        intersectPoint.x = ((side.b - ray.b) / (ray.a - side.a));
        intersectPoint.y = side.a * intersectPoint.x + side.b;
    } else if (ray.vertical && !side.vertical) {
        intersectPoint.x = ray.v[0].x;
        intersectPoint.y = side.a * intersectPoint.x + side.b;
    } else if (!ray.vertical && side.vertical) {
        intersectPoint.x = side.v[0].x;
        intersectPoint.y = ray.a * intersectPoint.x + ray.b;
    } else {
        return false;
    }
    if (is_inside(side, intersectPoint) && is_inside(ray, intersectPoint)) {
        return true;
    }
    return false;
}

bounding_box create_bounding_box() {
    bounding_box bb;
    bb.xMax = std::numeric_limits<double>::min();
    bb.xMin = std::numeric_limits<double>::max();
    bb.yMax = std::numeric_limits<double>::min();
    bb.yMin = std::numeric_limits<double>::max();
    return bb;
}

void update_bounding_box(Vector3 point, bounding_box *bb) {
    /* Setting the bounding box: */
    if (point.x > bb->xMax) {
        bb->xMax = point.x;
    } else if (point.x < bb->xMin) {
        bb->xMin = point.x;
    }
    if (point.y > bb->yMax) {
        bb->yMax = point.y;
    } else if (point.y < bb->yMin) {
        bb->yMin = point.y;
    }
}

bool insided_bounding_box(Vector3 point, bounding_box bb) {
    if ((point.x < bb.xMin) || (point.x > bb.xMax) || (point.y < bb.yMin) || (point.y > bb.yMax)) {
        return false;
    }
    return true;
}

bool contains(Vector3 point, bounding_box bb, std::vector<LineSegment> sides, int index) {
    if (insided_bounding_box(point, bb)) {
        LineSegment ray = create_ray(point, bb, index);
        int intersection = 0;
        for (int i = 0; i < sides.size(); i++) {
            if ((sides.at(i).index != index) && ray_intersect(ray, sides.at(i))) {
                intersection++;
            }
        }
        /* If the number of intersections is odd, then the point is inside the polygon: */
        if ((intersection % 2) == 1) {
            return true;
        }
    }
    return false;
}

void add_point(Vector3 p1, Vector3 p2, std::vector<LineSegment> &t, bounding_box *bb, bool first, int index) {
    if (first) {
        update_bounding_box(p1, bb);
    }
    update_bounding_box(p2, bb);
    LineSegment line(p1, p2, index);
    t.push_back(line);
}

void ray_casting(std::vector<contour> &polygons) {

    std::vector<LineSegment> segments;

    bounding_box bb = create_bounding_box();

    /*Creating the line segments of each contour: */
    for (int i = 0; i < polygons.size(); i++) {
        double area = 0.0;
        std::vector<Vector3> Pi = polygons.at(i).P;
        for (int j = 1; j < Pi.size(); j++) {
            Vector3 p0 = Pi.at(j - 1);
            Vector3 p1 = Pi.at(j + 0);
            area += (p0.x * p1.y - p0.y * p1.x);
            add_point(p0, p1, segments, &bb, (j == 1 ? true : false), i);
            if (j == Pi.size() - 1) {
                add_point(p1, Pi.at(0), segments, &bb, (j == 1 ? true : false), i);
                area += (p1.x * Pi.at(0).y - p1.y * Pi.at(0).x);
            }
        }
        area /= 2.0;
        if (area < 0.0) {
            polygons.at(i).clockwise = true;
        } else {
            polygons.at(i).clockwise = false;
        }
    }

    /*Using the point in polygon algorithm to test the first segment of each contour: */
    for (int i = 0; i < polygons.size(); i++) {
        std::vector<Vector3> Pi = polygons.at(i).P;
        if (contains(Pi.at(0), bb, segments, i)) {
            /*Internal contour: */
            polygons.at(i).external = false;
        } else {
            /*External contour: */
            polygons.at(i).external = true;
        }

        /*Reversing contours: */
        if (polygons.at(i).external && polygons.at(i).clockwise) {
            std::reverse(polygons.at(i).P.begin(), polygons.at(i).P.end());
            polygons.at(i).clockwise = false;
        } else if (!polygons.at(i).external && !polygons.at(i).clockwise) {
            std::reverse(polygons.at(i).P.begin(), polygons.at(i).P.end());
            polygons.at(i).clockwise = true;
        }
    }
    segments.clear();
}

void TrivialSlicing(const TriangleMesh &mesh, std::vector<float> &planes) {

    int numberPlanes = planes.size(); /* Number of planes. */
    const std::vector<Triangle> &T = mesh.getTriangles();
    std::vector<LineSegment> segs[numberPlanes];

    /* Enumerate the slicing planes: */
    for (int p = 0; p < numberPlanes; p++) {
        /* Enumerate all triangles of the mesh:*/
        for (auto it = T.begin(), itEnd = T.end(); it != itEnd; ++it) {
            Triangle t = *it; /*Current triangle.*/
            /*Does the plane intersect the triangle?:*/
            if ((t.zMin < planes[p]) && (t.zMax > planes[p])) {
                /* Compute and save the intersection: */
                LineSegment seg = R3_Mesh_Triangle_slice(t, planes[p]);
                seg.v[0].x = round(seg.v[0].x * 100.0) / 100.0;
                seg.v[0].y = round(seg.v[0].y * 100.0) / 100.0;
                seg.v[1].x = round(seg.v[1].x * 100.0) / 100.0;
                seg.v[1].y = round(seg.v[1].y * 100.0) / 100.0;
                if (seg.v[0].distTo(seg.v[1]) > 0.0001) {
                    segs[p].push_back(seg);
                }
            }
        }
        if (!segs[p].empty()) {
            std::vector<contour> contours;
            TrivialLoopClosure(segs[p], contours);
            ray_casting(contours);
            segs[p].clear();
//            cout << planes[p] << " " << contours.size() << endl;
//            for (auto contour: contours) {
//                cout << "Clockwise: " << contour.clockwise << endl;
//                for (auto p: contour.P) {
//                    cout << p << " " << endl;
//                }
//                cout << endl;
//            }
        }
    }
}

int main() {
    TriangleMesh mesh;
    double thickness{0.1};

    stlToMeshInMemory("../cube_10x10x10.stl", mesh);
//    stlToMeshInMemory("../pisa.stl", mesh);
//    stlToMeshInMemory("../bblocky.stl", mesh);
    std::vector<float> planes = compute_planes(mesh, thickness);
    std::cout << "Planes: " << planes.size() << " [";
    for (auto z: planes) {
        std::cout << z << " ";
    }
    std::cout << std::endl;

    Timer timer;

    timer.reset();
    TrivialSlicing(mesh, planes);
    std::cout << "Slicing: " << timer.elapsed() << " seconds" << std::endl;

    return 0;
}