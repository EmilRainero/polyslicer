#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include "Timer.h"


using namespace std;


double eps = 0.0001; // round all vertices in stl file to wthin this epsilon

class v3 {

public:

    v3 (float _x=0, float _y=0, float _z=0) : x(_x), y(_y), z(_z) {}

    float distTo (const v3 &pt) const {
        return sqrt ( pow(fabs(x-pt.x), 2.0) + pow(fabs(y-pt.y), 2.0) + pow(fabs(z-pt.z), 2.0) );
    }

//    array<float,3> getCoords() {
//        array<float,3> c = {{x, y, z}};
//        return c;
//    }

    float dotproduct (const v3 &v) const {
        return (x*v.x + y*v.y + z*v.z);
    }

//    void transform (const glm::mat4 &mat) {
//        glm::vec4 v = glm::vec4(x, y, z, 1.0);
//        glm::vec4 vt = mat*v;
//        x = (vt.x); y = (vt.y); z = (vt.z);
//    }

    v3& operator-=(const v3 &pt) {
        x = (x-pt.x);
        y = (y-pt.y);
        z = (z-pt.z);
        return *this;
    }

    v3 operator-(const v3 &pt) {
        return v3 ((x-pt.x), (y-pt.y), (z-pt.z));
    }

    v3 operator+(const v3 &pt) {
        return v3 ((x+pt.x), (y+pt.y), (z+pt.z));
    }

    v3 operator/(float a) {
        return v3 ((x/a), (y/a), (z/a));
    }

    v3 operator*(float a) {
        return v3 ((x*a), (y*a), (z*a));
    }

    bool operator<(const v3 &pt) const {
        return z < pt.z;
    }

    bool operator>(const v3 &pt) const {
        return z > pt.z;
    }

    bool operator==(const v3 &pt) const {
        return distTo(pt) < 0.005;
    }

    bool operator!=(const v3 &pt) const {
        return distTo(pt) > 0.005;
    }

    float normalize() const {
        return sqrt(x*x+y*y+z*z);
    }

    string getLabel() const {
        stringstream ss;
        ss << x << "|" << y << "|" << z;
        return ss.str();
    }

    friend ostream& operator<<(ostream& os, const v3& v) {
        os << "x: " << v.x << "; y: " << v.y << "; z: " << v.z;
        return os;
    }

public:

    float x;
    float y;
    float z;
};

class Triangle    {

public:

    Triangle(v3 n, v3 v0, v3 v1, v3 v2) : normal(n) {
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        zMin = +99999999.9;
        zMax = -99999999.9;
        setZMin(v0.z); setZMin(v1.z); setZMin(v2.z);
        setZMax(v0.z); setZMax(v1.z); setZMax(v2.z);
    }

    void setZMin (float z) {
        if (z < zMin) {
            zMin = z;
        }
    }

    void setZMax (float z) {
        if (z > zMax) {
            zMax = z;
        }
    }

    Triangle& operator-=(const v3 &pt) {
        v[0] -= pt;
        v[1] -= pt;
        v[2] -= pt;
        return *this;
    }

    bool operator<(const Triangle &t) {
        return zMin < t.zMin;
    }

    friend ostream& operator<<(ostream& os, const Triangle& t) {
        os << "V0: (" << t.v[0] << "); V1: (" << t.v[1] << "); V2: (" << t.v[2] << ")";
        return os;
    }

    bool degenerate () {
        if (v[0].distTo(v[1]) < 0.000001) { return true; }
        if (v[1].distTo(v[2]) < 0.000001) { return true; }
        if (v[2].distTo(v[0]) < 0.000001) { return true; }
        return false;
    }


public:

    v3 v[3];
    v3 normal;
    float zMin;
    float zMax;
};

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

    v3 meshAABBSize() const {
        return v3 ( topRightVertex.x - bottomLeftVertex.x,
                    topRightVertex.y - bottomLeftVertex.y,
                    topRightVertex.z - bottomLeftVertex.z );
    }

    const vector<Triangle>& getTriangles() const { return triangles; }

    v3 getBottomLeftVertex() const { return bottomLeftVertex; }

    v3 getTopRightVertex() const { return topRightVertex; }

public:

    int meshSize;
    v3 bottomLeftVertex;
    v3 topRightVertex;
    vector<Triangle> triangles;
};

class LineSegment {

public:

    LineSegment (v3 p0=v3(), v3 p1=v3(), int i=0) {
        v[0] = p0;
        v[1] = p1;
        index = i;
        vertical = false;
        if ((v[1].x - v[0].x) != 0) {
            a = (v[1].y - v[0].y)/(v[1].x - v[0].x);
            b = (v[0].y - (a * v[0].x));
        }
        else {
            vertical = true;
        }
    }

    bool operator==(const LineSegment &ls) const {
        return ((v[0] == ls.v[0]) && (v[1] == ls.v[1]));
    }

    friend ostream& operator<<(ostream& os, const LineSegment& ls) {
        os << "V0: (" << ls.v[0] << "); V1: (" << ls.v[1] << ")";
        return os;
    }

public:

    v3 v[2];
    double a;
    double b;
    bool vertical;
    int index;
};



float xround (float x, double eps, int mod, int rem) {
    double y = round((double)x/(mod * eps));
    double z = (y * mod + rem) * eps;
    return (float)z;
}

v3 v3_round (float x, float y, float z, double eps) {
    v3 p;
    p.x = xround(x, eps, 2, 0);
    p.y = xround(y, eps, 2, 0);
    p.z = xround(z, eps, 2, 0);
    return p;
}

Triangle make_triangle (
        float n0, float n1, float n2,
        float f0, float f1, float f2,
        float f3, float f4, float f5,
        float f6, float f7, float f8,
        double eps)
{

    return Triangle (v3(n0, n1, n2), v3_round(f0, f1, f2, eps), v3_round(f3, f4, f5, eps), v3_round(f6, f7, f8, eps));
}

int stlToMeshInMemory (const char *stlFile, TriangleMesh& mesh) {

    int ndegenerated = 0;

    FILE *f = fopen (stlFile, "rb");
    if (!f) {
        return 1;
    }
    char title[80];
    int nFaces;
    int err;
    err = fread (title, 80, 1, f);
    err = fread ((void*)&nFaces, 4, 1, f);
    float v[12]; /* normal = 3, vertices = 9 (12) */
    unsigned short uint16;
    /* Every Face is 50 Bytes: Normal(3*float), Vertices(9*float), 2 Bytes Spacer */
    for (size_t i=0; i<nFaces; ++i) {
        for (size_t j=0; j<12; ++j) {
            err = fread((void*)&v[j], sizeof(float), 1, f);
        }
        err = fread((void*)&uint16, sizeof(unsigned short), 1, f); // spacer between successive faces
        Triangle t = make_triangle (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], eps);
        if (!t.degenerate()) {
            mesh.push_back (t);
        }
        else {
            ndegenerated++;
        }
    }
    fclose(f);

    cout << "Number triangles: " << nFaces << endl;
    cout << "Number degenerated triangles: " << ndegenerated << endl;
    return 0;
}

vector<float> compute_planes(TriangleMesh& mesh, double thickness) {
    vector<float> planes;

    double model_zmax = mesh.getTopRightVertex().z;
    double model_zmin = mesh.getBottomLeftVertex().z;

    double rounded_thickness = xround (thickness, eps, 2, 0); /*Plane spacing even multiple of {eps}*/
    double P0 = xround (model_zmin - rounded_thickness, eps, 2, 1); /*First plane odd multiple of {eps}.*/

    int no_planes = 1 + (int)((model_zmax + rounded_thickness - P0)/rounded_thickness); // includes padding

    cout << "Model Z: [" << model_zmin << "  " << model_zmax << "] = " << model_zmax - model_zmin << " height" << endl;

    for (size_t i = 0; i < no_planes; i++) {
        /* Building the vector with the slice z coordinates: */
        float Pi = (float)(P0 + i * rounded_thickness);
        if ((Pi > model_zmin) && (Pi < model_zmax)) {
            planes.push_back (Pi);
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

v3 R3_Mesh_Side_slice (v3 vi, v3 vj, float Z) {
    double dx = vj.x - vi.x;
    double dy = vj.y - vi.y;
    double dz = vj.z - vi.z;
//    assert(dz != 0);
    double frac = (Z - vi.z)/dz;
    float xint = (float)(frac*dx + (double)vi.x);
    float yint = (float)(frac*dy + (double)vi.y);
    return (v3){ .x = xint, .y = yint, .z = Z };
}

LineSegment R3_Mesh_Triangle_slice (Mesh_Triangle_Node_t *t, float Z) {
//    assert((t->t.zMin < Z) && (t->t.zMax > Z));
    int np = 0; /* Number of segment endpoints found */
    LineSegment seg;
    for (int i = 0; i < 3; i++) {
        /* Get side {i} of triangle: */
        int j = (i == 2 ? 0 : i+1);
        v3 vi = (t->t.v[i]);
        v3 vj = (t->t.v[j]);
        /* Check for intersection of plane with {vi--vj}. */
        /* Must consider segment closed at bottom and open at top in case {Z} goes through a vertex. */
        float vzMin = (vi.z < vj.z ? vi.z : vj.z);
        float vzMax = (vi.z > vj.z ? vi.z : vj.z);
        if ((vzMin <= Z) && (vzMax > Z)) {
            v3 p = R3_Mesh_Side_slice (vi, vj, Z);
//            assert(np < 2);
            seg.v[np] = p;
            np++;
        }
    }
//    assert(np == 2);
    return seg;
}

LineSegment R3_Mesh_Triangle_slice (Triangle t, float Z) {
//    assert((t.zMin < Z) && (t.zMax > Z));
    int np = 0; /* Number of segment endpoints found */
    LineSegment seg;
    for (int i = 0; i < 3; i++) {
        /* Get side {i} of triangle: */
        int j = (i == 2 ? 0 : i+1);
        v3 vi = (t.v[i]);
        v3 vj = (t.v[j]);
        /* Check for intersection of plane with {vi--vj}. */
        /* Must consider segment closed at bottom and open at top in case {Z} goes through a vertex. */
        float vzMin = (vi.z < vj.z ? vi.z : vj.z);
        float vzMax = (vi.z > vj.z ? vi.z : vj.z);
        if ((vzMin <= Z) && (vzMax > Z)) {
            v3 p = R3_Mesh_Side_slice (vi, vj, Z);
//            assert(np < 2);
            seg.v[np] = p;
            np++;
        }
    }
//    assert(np == 2);
    return seg;
}

v3 TrivialClosure_find_A (vector<LineSegment>& segs, v3 u) {
    for (vector<LineSegment>::iterator j = segs.begin(); j != segs.end(); ++j) {
        v3 p0 = (*j).v[0];
        v3 p1 = (*j).v[1];
        if ((p0.x == u.x) && (p0.y == u.y)) {
            segs.erase(j);
            return p1;
        } else if ((p1.x == u.x) && (p1.y == u.y)) {
            segs.erase(j);
            return p0;
        }
    }
    return {+INFINITY,+INFINITY,+INFINITY};
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
    vector<v3> P;
} contour;

void TrivialLoopClosure (vector<LineSegment> segs, vector<contour>& contours) {

    while (segs.size() > 0) {

        /* Get another contour: */
        vector<v3> P;

        /* Get the first segment: */
        v3 last;
        v3 current;
        v3 prev;
        {
            std::vector<LineSegment>::iterator i = segs.begin();
            last = (*i).v[0];
            current = (*i).v[1];
            {
                P.push_back (last);
            }
            prev = last;
            segs.erase(i);
        }

        /* Get additional segments until loop is closed: */
        bool open = false;
        do {
            /* Find and delete another segment with endpoint {current}, advance {prev,current} */
            v3 next = TrivialClosure_find_A (segs, current);
            //v3 next = TrivialClosure_find_B (segs, current, plane);
            if ((next.x == +INFINITY) || (next.y == +INFINITY)) {
                /* Open contour?! */
                open = true;
                break;
            }
            {
                P.push_back (current);
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

LineSegment create_ray (v3 point, bounding_box bb, int index) {
    /* Create outside point: */
    double epsilon = (bb.xMax - bb.xMin) / 100.0;
    v3 outside (bb.xMin - epsilon, bb.yMin);
    LineSegment v (outside, point, index);
    return v;
}

bool is_inside (LineSegment line, v3 point) {
    double maxX = (line.v[0].x > line.v[1].x) ? line.v[0].x : line.v[1].x;
    double minX = (line.v[0].x < line.v[1].x) ? line.v[0].x : line.v[1].x;
    double maxY = (line.v[0].y > line.v[1].y) ? line.v[0].y : line.v[1].y;
    double minY = (line.v[0].y < line.v[1].y) ? line.v[0].y : line.v[1].y;
    if ((point.x >= minX && point.x <= maxX) && (point.y >= minY && point.y <= maxY)) {
        return true;
    }
    return false;
}

bool ray_intersect (LineSegment ray, LineSegment side) {
    v3 intersectPoint;
    /* If both vectors aren't from the kind of x=1 lines then go into: */
    if (!ray.vertical && !side.vertical) {
        /* Check if both vectors are parallel. If they are parallel then no intersection point will exist: */
        if (ray.a - side.a == 0) {
            return false;
        }
        intersectPoint.x = ((side.b - ray.b) / (ray.a - side.a));
        intersectPoint.y = side.a * intersectPoint.x + side.b;
    }
    else if (ray.vertical && !side.vertical) {
        intersectPoint.x = ray.v[0].x;
        intersectPoint.y = side.a * intersectPoint.x + side.b;
    }
    else if (!ray.vertical && side.vertical) {
        intersectPoint.x = side.v[0].x;
        intersectPoint.y = ray.a * intersectPoint.x + ray.b;
    }
    else {
        return false;
    }
    if (is_inside(side, intersectPoint) && is_inside(ray, intersectPoint)) {
        return true;
    }
    return false;
}

bounding_box create_bounding_box () {
    bounding_box bb;
    bb.xMax = std::numeric_limits<double>::min();
    bb.xMin = std::numeric_limits<double>::max();
    bb.yMax = std::numeric_limits<double>::min();
    bb.yMin = std::numeric_limits<double>::max();
    return bb;
}

void update_bounding_box (v3 point, bounding_box *bb) {
    /* Setting the bounding box: */
    if (point.x > bb->xMax) {
        bb->xMax = point.x;
    }
    else if (point.x < bb->xMin) {
        bb->xMin = point.x;
    }
    if (point.y > bb->yMax) {
        bb->yMax = point.y;
    }
    else if (point.y < bb->yMin) {
        bb->yMin = point.y;
    }
}

bool insided_bounding_box (v3 point, bounding_box bb) {
    if ( (point.x < bb.xMin) || (point.x > bb.xMax) || (point.y < bb.yMin) || (point.y > bb.yMax) ) {
        return false;
    }
    return true;
}

bool contains (v3 point, bounding_box bb, vector<LineSegment> sides, int index) {
    if (insided_bounding_box(point, bb)) {
        LineSegment ray = create_ray (point, bb, index);
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
void add_point (v3 p1, v3 p2, vector<LineSegment> &t, bounding_box *bb, bool first, int index) {
    if (first) {
        update_bounding_box (p1, bb);
    }
    update_bounding_box (p2, bb);
    LineSegment line (p1, p2, index);
    t.push_back(line);
}

void ray_casting (vector<contour> &polygons) {

    vector<LineSegment> segments;

    bounding_box bb = create_bounding_box ();

    /*Creating the line segments of each contour: */
    for (int i = 0; i < polygons.size(); i++) {
        double area = 0.0;
        vector<v3> Pi = polygons.at(i).P;
        for (int j = 1; j < Pi.size(); j++) {
            v3 p0 = Pi.at(j-1);
            v3 p1 = Pi.at(j+0);
            area += (p0.x * p1.y - p0.y * p1.x);
            add_point (p0, p1, segments, &bb, (j == 1 ? true : false), i);
            if (j == Pi.size()-1) {
                add_point (p1, Pi.at(0), segments, &bb, (j == 1 ? true : false), i);
                area += (p1.x * Pi.at(0).y - p1.y * Pi.at(0).x);
            }
        }
        area /= 2.0;
        if (area < 0.0) {
            polygons.at(i).clockwise = true;
        }
        else {
            polygons.at(i).clockwise = false;
        }
    }

    /*Using the point in polygon algorithm to test the first segment of each contour: */
    for (int i = 0; i < polygons.size(); i++) {
        vector<v3> Pi = polygons.at(i).P;
        if (contains (Pi.at(0), bb, segments, i)) {
            /*Internal contour: */
            polygons.at(i).external = false;
        }
        else {
            /*External contour: */
            polygons.at(i).external = true;
        }

        /*Reversing contours: */
        if (polygons.at(i).external && polygons.at(i).clockwise) {
            std::reverse(polygons.at(i).P.begin(), polygons.at(i).P.end());
            polygons.at(i).clockwise = false;
        }
        else if (!polygons.at(i).external && !polygons.at(i).clockwise) {
            std::reverse(polygons.at(i).P.begin(), polygons.at(i).P.end());
            polygons.at(i).clockwise = true;
        }
    }
    segments.clear();
}

void TrivialSlicing (const TriangleMesh& mesh, vector<float>& planes) {

    int numberPlanes = planes.size(); /* Number of planes. */

    const vector<Triangle> &T = mesh.getTriangles();

    vector<LineSegment> segs[numberPlanes];

    int intersections = 0;

    /* Enumerate the slicing planes: */
    for (int p = 0; p < numberPlanes; p++) {
        /* Enumerate all triangles of the mesh:*/
        int intersections_per_plane = 0;
        for (auto it = T.begin(), itEnd = T.end(); it != itEnd; ++it) {
            Triangle t = *it; /*Current triangle.*/
            /*Does the plane intersect the triangle?:*/
            if ((t.zMin < planes[p]) && (t.zMax > planes[p])) {
                /* Compute and save the intersection: */
                LineSegment seg = R3_Mesh_Triangle_slice (t, planes[p]);
                seg.v[0].x = round(seg.v[0].x * 100.0) / 100.0;
                seg.v[0].y = round(seg.v[0].y * 100.0) / 100.0;
                seg.v[1].x = round(seg.v[1].x * 100.0) / 100.0;
                seg.v[1].y = round(seg.v[1].y * 100.0) / 100.0;
                if (seg.v[0].distTo(seg.v[1]) > 0.0001) {
                    segs[p].push_back(seg);
                }
                intersections++;
                intersections_per_plane++;
            }
        }
//        cout << planes[p] << " " << segs[p].size() << "  " << intersections_per_plane << endl;
        if (!segs[p].empty()) {
            vector<contour> contours;
            TrivialLoopClosure (segs[p], contours);
            ray_casting (contours);
            segs[p].clear();
//            cout << planes[p] << " " << contours.size() << endl;
            for (auto contour: contours) {
//                cout << "Clockwise: " << contour.clockwise << endl;
                for (auto p: contour.P) {
//                    cout << p << " " << endl;
                }
//                cout << endl;
            }
        }
    }
    cout << "Intersections: " << intersections << endl;

//    if (chaining) {
//        /*Loop-Closure:*/
//        clock_t contour_begin = clock();
//        for (size_t p = 0; p < numberPlanes; p++) {
//
//        }
//        clock_t contour_end = clock();
//        loopclosure_time = double(contour_end - contour_begin)/CLOCKS_PER_SEC;
//        /*End-Loop-Closure*/
//    }
//    else {
//        export_svg_no_chaining ("segments.svg", segs, k, mesh->meshAABBSize());
//    }
}

int main() {
    TriangleMesh mesh;
    double thickness{0.1};

//    stlToMeshInMemory("../cube_10x10x10.stl", mesh);
    stlToMeshInMemory("../pisa.stl", mesh);
//    stlToMeshInMemory("../bblocky.stl", mesh);
    vector<float> planes = compute_planes(mesh, thickness);
    cout << "Planes: " << planes.size() << " [";
    for (auto z: planes) {
        cout << z << " ";
    }
    cout << endl;

    Timer timer;

    timer.reset();
    TrivialSlicing(mesh, planes);
    cout << "Slicing: " << timer.elapsed() << " seconds" << endl;

    return 0;
}