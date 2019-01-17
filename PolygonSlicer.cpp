#include <limits>
#include <algorithm>

#include "PolygonSlicer.h"
#include "LineSegment.h"
#include "Contour.h"
#include "BoundingBox.h"
#include "Timer.h"


float PolygonSlicer::roundToNearestEpsilon(float x, double epsilon, int mod, int rem) {
    double y = round((double) x / (mod * epsilon));
    double z = (y * mod + rem) * epsilon;
    return (float) z;
}

Vector3 PolygonSlicer::roundVector3(float x, float y, float z, double epsilon) {
    Vector3 p;
    p.x = roundToNearestEpsilon(x, epsilon, 2, 0);
    p.y = roundToNearestEpsilon(y, epsilon, 2, 0);
    p.z = roundToNearestEpsilon(z, epsilon, 2, 0);
    return p;
}

Triangle PolygonSlicer::makeRoundedTriangle(
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

std::vector<float> PolygonSlicer::computePlanes(TriangleMesh &mesh, double thickness, double epsilon) {
    std::vector<float> planes;

    double model_zmax = mesh.getTopRightVertex().z;
    double model_zmin = mesh.getBottomLeftVertex().z;

    double rounded_thickness = roundToNearestEpsilon(thickness, epsilon, 2, 0); /*Plane spacing even multiple of {epsilon}*/
    double P0 = roundToNearestEpsilon(model_zmin - rounded_thickness, epsilon, 2, 1); /*First plane odd multiple of {epsilon}.*/

    int no_planes = 1 + (int) ((model_zmax + rounded_thickness - P0) / rounded_thickness);

//    std::cout << "Model Z: [" << model_zmin << "  " << model_zmax << "] = " << model_zmax - model_zmin << " height"
//              << std::endl;

    for (size_t i = 0; i < no_planes; i++) {
        float Pi = (float) (P0 + i * rounded_thickness);
        if ((Pi > model_zmin) && (Pi < model_zmax)) {
            planes.push_back(Pi);
        }
    }

    return planes;
}


PolygonSlicer::Mesh_Triangle_List_t* PolygonSlicer::Mesh_Triangle_List_create (void) {
    Mesh_Triangle_List_t *L = (Mesh_Triangle_List_t *)malloc(sizeof(Mesh_Triangle_List_t));
    L->head = NULL;
    L->tail = NULL;
    return L;
}

void PolygonSlicer::Mesh_Triangle_List_insert (Triangle t, Mesh_Triangle_List_t *L) {
    Mesh_Triangle_Node_t *node = (Mesh_Triangle_Node_t *)malloc(sizeof(Mesh_Triangle_Node_t));
    node->t = t;
    node->next = L->head;
    node->prev = NULL;
    if (L->head == NULL) {
        /*New head*/
        L->head = L->tail = node;
    }
    else {
        L->head->prev = node;
        L->head = node;
    }
}

/*-----------------------------------------------------------------------*/
void PolygonSlicer::Mesh_Triangle_List_union (Mesh_Triangle_List_t *L1, Mesh_Triangle_List_t *L2) {
    if ( (L1->head != NULL) && (L2->head != NULL) ) {
        L1->tail->next = L2->head;
        L2->head->prev = L1->tail;
        L1->tail = L2->tail;;
    }
    else if (L2->head != NULL) {
        L1->head = L2->head;
        L1->tail = L2->tail;
    }
}

/*-----------------------------------------------------------------------*/
PolygonSlicer::Mesh_Triangle_Node_t* PolygonSlicer::Mesh_Triangle_List_remove (Mesh_Triangle_List_t *L, Mesh_Triangle_Node_t *node) {
    if ((node->prev == NULL) && (node->next == NULL)) {
        free (node);
        L->head = NULL;
        L->tail = NULL;
        return NULL;
    }
    else if (node->prev == NULL) {
        node->next->prev = NULL;
        L->head = node->next;
        free (node);
        return L->head;
    }
    else if (node->next == NULL) {
        node->prev->next = NULL;
        L->tail = node->prev;
        free (node);
        return NULL;
    }
    else {
        Mesh_Triangle_Node_t *next = node->next;
        node->next->prev = node->prev;
        node->prev->next = next;
        free (node);
        return next;
    }
}

Vector3 PolygonSlicer::R3_Mesh_Side_slice(Vector3 vi, Vector3 vj, float Z) {
    double dx = vj.x - vi.x;
    double dy = vj.y - vi.y;
    double dz = vj.z - vi.z;
//    assert(dz != 0);
    double frac = (Z - vi.z) / dz;
    float xint = (float) (frac * dx + (double) vi.x);
    float yint = (float) (frac * dy + (double) vi.y);
    Vector3 result(xint, yint, Z);
    return result;
}


LineSegment PolygonSlicer::R3_Mesh_Triangle_slice (Mesh_Triangle_Node_t *t, float Z) {
//    assert((t->t.zMin < Z) && (t->t.zMax > Z));
    int np = 0; /* Number of segment endpoints found */
    LineSegment seg;
    for (int i = 0; i < 3; i++) {
        /* Get side {i} of triangle: */
        int j = (i == 2 ? 0 : i+1);
        Vector3 vi = (t->t.v[i]);
        Vector3 vj = (t->t.v[j]);
        /* Check for intersection of plane with {vi--vj}. */
        /* Must consider segment closed at bottom and open at top in case {Z} goes through a vertex. */
        float vzMin = (vi.z < vj.z ? vi.z : vj.z);
        float vzMax = (vi.z > vj.z ? vi.z : vj.z);
        if ((vzMin <= Z) && (vzMax > Z)) {
            Vector3 p = R3_Mesh_Side_slice (vi, vj, Z);
//            assert(np < 2);
            seg.v[np] = p;
            np++;
        }
    }
//    assert(np == 2);
    return seg;
}

LineSegment PolygonSlicer::R3_Mesh_Triangle_slice(Triangle t, float Z) {
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

Vector3 PolygonSlicer::TrivialClosure_find_A(std::vector<LineSegment> &segs, Vector3 u) {
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

void PolygonSlicer::TrivialLoopClosure(std::vector<LineSegment> lineSegments, std::vector<Contour>& contour) {

    while (lineSegments.size() > 0) {

        /* Get another contour: */
        std::vector<Vector3> P;

        Vector3 last;
        Vector3 current;
        Vector3 prev;
        {
            std::vector<LineSegment>::iterator i = lineSegments.begin();
            last = (*i).v[0];
            current = (*i).v[1];
            {
                P.push_back(last);
            }
            prev = last;
            lineSegments.erase(i);
        }

        /* Get additional segments until loop is closed: */
        bool open = false;
        do {
            /* Find and delete another segment with endpoint {current}, advance {prev,current} */
            Vector3 next = TrivialClosure_find_A(lineSegments, current);
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
        contour.push_back({false, false, P});
    }
}

LineSegment PolygonSlicer::create_ray(Vector3 point, BoundingBox bbox, int index) {
    /* Create outside point: */
    double epsilon = (bbox.xMax - bbox.xMin) / 100.0;
    Vector3 outside(bbox.xMin - epsilon, bbox.yMin);
    LineSegment v(outside, point, index);
    return v;
}

bool PolygonSlicer::is_inside(LineSegment line, Vector3 point) {
    double maxX = (line.v[0].x > line.v[1].x) ? line.v[0].x : line.v[1].x;
    double minX = (line.v[0].x < line.v[1].x) ? line.v[0].x : line.v[1].x;
    double maxY = (line.v[0].y > line.v[1].y) ? line.v[0].y : line.v[1].y;
    double minY = (line.v[0].y < line.v[1].y) ? line.v[0].y : line.v[1].y;
    if ((point.x >= minX && point.x <= maxX) && (point.y >= minY && point.y <= maxY)) {
        return true;
    }
    return false;
}

bool PolygonSlicer::ray_intersect(LineSegment ray, LineSegment side) {
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

void PolygonSlicer::update_bounding_box(Vector3 point, BoundingBox *bbox) {
    /* Setting the bounding box: */
    if (point.x > bbox->xMax) {
        bbox->xMax = point.x;
    } else if (point.x < bbox->xMin) {
        bbox->xMin = point.x;
    }
    if (point.y > bbox->yMax) {
        bbox->yMax = point.y;
    } else if (point.y < bbox->yMin) {
        bbox->yMin = point.y;
    }
}

bool PolygonSlicer::insided_bounding_box(Vector3 point, BoundingBox bb) {
    if ((point.x < bb.xMin) || (point.x > bb.xMax) || (point.y < bb.yMin) || (point.y > bb.yMax)) {
        return false;
    }
    return true;
}

bool PolygonSlicer::contains(Vector3 point, BoundingBox bb, std::vector<LineSegment> sides, int index) {
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

void PolygonSlicer::add_point(Vector3 p1, Vector3 p2, std::vector<LineSegment> &t, BoundingBox *bbox, bool first, int index) {
    if (first) {
        update_bounding_box(p1, bbox);
    }
    update_bounding_box(p2, bbox);
    LineSegment line(p1, p2, index);
    t.push_back(line);
}

void PolygonSlicer::ray_casting(std::vector<Contour> &polygons) {

    std::vector<LineSegment> segments;

    BoundingBox bbox;

    /*Creating the line segments of each contour: */
    for (int i = 0; i < polygons.size(); i++) {
        double area = 0.0;
        std::vector<Vector3> Pi = polygons.at(i).points;
        for (int j = 1; j < Pi.size(); j++) {
            Vector3 p0 = Pi.at(j - 1);
            Vector3 p1 = Pi.at(j + 0);
            area += (p0.x * p1.y - p0.y * p1.x);
            add_point(p0, p1, segments, &bbox, (j == 1 ? true : false), i);
            if (j == Pi.size() - 1) {
                add_point(p1, Pi.at(0), segments, &bbox, (j == 1 ? true : false), i);
                area += (p1.x * Pi.at(0).y - p1.y * Pi.at(0).x);
            }
        }
        area /= 2.0;
        if (area >= 0.0) {
            polygons.at(i).clockwise = true;
        } else {
            polygons.at(i).clockwise = false;
        }
    }

    /*Using the point in polygon algorithm to test the first segment of each contour: */
    for (int i = 0; i < polygons.size(); i++) {
        std::vector<Vector3> Pi = polygons.at(i).points;
        if (contains(Pi.at(0), bbox, segments, i)) {
            /*Internal contour: */
            polygons.at(i).external = false;
        } else {
            /*External contour: */
            polygons.at(i).external = true;
        }

        /*Reversing contours: */
        if (polygons.at(i).external && polygons.at(i).clockwise) {
            std::reverse(polygons.at(i).points.begin(), polygons.at(i).points.end());
            polygons.at(i).clockwise = false;
        } else if (!polygons.at(i).external && !polygons.at(i).clockwise) {
            std::reverse(polygons.at(i).points.begin(), polygons.at(i).points.end());
            polygons.at(i).clockwise = true;
        }
    }
    segments.clear();
}

std::vector<Layer *> PolygonSlicer::TrivialSlicing(const TriangleMesh &mesh, std::vector<float> &planes) {

    int numberPlanes = planes.size();
    const std::vector<Triangle> &triangles = mesh.getTriangles();
    std::vector<LineSegment> lineSegments[numberPlanes];
    std::vector<Layer *> layers;

    for (int p = 0; p < numberPlanes; p++) {
        float z = planes[p];
        layers.push_back(new Layer(z));

        for (auto it = triangles.begin(), itEnd = triangles.end(); it != itEnd; ++it) {
            Triangle t = *it; /*Current triangle.*/
            /*Does the plane intersect the triangle?:*/
            if ((t.zMin < z) && (t.zMax > z)) {
                /* Compute and save the intersection: */
                LineSegment lineSegment = R3_Mesh_Triangle_slice(t, z);
                lineSegment.v[0].x = round(lineSegment.v[0].x * 100.0) / 100.0;
                lineSegment.v[0].y = round(lineSegment.v[0].y * 100.0) / 100.0;
                lineSegment.v[1].x = round(lineSegment.v[1].x * 100.0) / 100.0;
                lineSegment.v[1].y = round(lineSegment.v[1].y * 100.0) / 100.0;
                if (lineSegment.v[0].distTo(lineSegment.v[1]) > 0.0001) {
                    lineSegments[p].push_back(lineSegment);
                }
            }
        }
        if (!lineSegments[p].empty()) {
            TrivialLoopClosure(lineSegments[p], layers[p]->contours);
            ray_casting(layers[p]->contours);
            lineSegments[p].clear();
        }
    }

    return layers;
}

/* Assumes that {P[0..k-1]} is a list of {k} strictly increasing {Z}
  coordinates. Returns a vector of {k+1} lists {L[0..k]} such that {L[p]}
  contains all triangles of the {mesh} that have {zMin} between {P[p-1]}
  and {P[p]}, assuming that {P[-1] = -oo} and {P[k] = +oo}. If {delta > 0},
  assumes that {P[p]-P[p-1] = delta} for {p} in {1..k-1}. If {srt} is true,
  assumes that the triangles are already sorted by increasing {zMin}. */
PolygonSlicer::Mesh_Triangle_List_t** PolygonSlicer::IncrementalSlicing_buildLists (const TriangleMesh& mesh, std::vector<float> P, double delta) {

    int k = P.size(); /* Number of planes. */

    Mesh_Triangle_List_t **L = (Mesh_Triangle_List_t **)malloc((k+1) * sizeof(Mesh_Triangle_List_t *));

    for (size_t p = 0; p <= k; p++) { L[p] = Mesh_Triangle_List_create(); }

    const std::vector<Triangle> &T = mesh.triangles;

    int n = T.size(); /* Number of triangles. */

    /* Uniform slicing - compute list index: */
    for (auto it = T.begin(), itEnd = T.end(); it != itEnd; ++it) {
        Triangle t = *it;
        int p;
        if (t.zMin < P[0]) {
            p = 0;
        }
        else if (t.zMin > P[k-1]) {
            p = k;
        }
        else {
            p = floor((t.zMin - P[0])/delta) + 1;
        }
        Mesh_Triangle_List_insert (t, L[p]);
    }

    return L;
}

std::vector<Layer *> PolygonSlicer::IncrementalSlicing(const TriangleMesh& mesh, std::vector<float> &planes, double delta) {

    /*Slicing*/
    clock_t slice_begin = clock();

    int k = planes.size();

    std::vector <LineSegment> segs[k];

    /* Classify triangles by the plane gaps that contain their {zMin}: */
    Mesh_Triangle_List_t **L = IncrementalSlicing_buildLists(mesh, planes, delta);
    /* Now perform a plane sweep from bottom to top: */

    Mesh_Triangle_List_t *A = Mesh_Triangle_List_create(); /* Active triangle list. */
    for (int p = 0; p < k; p++) {
        /* Add triangles that start between {P[p-1]} and {P[p]}: */
        Mesh_Triangle_List_union(A, L[p]);
        /* Scan the active triangles: */
        Mesh_Triangle_Node_t *aux = A->head;
        while (aux != NULL) {
            Mesh_Triangle_Node_t *next = aux->next;
            if (aux->t.zMax < planes[p]) {
                /* Triangle is exhausted: */
                Mesh_Triangle_List_remove(A, aux);
            } else {
                /* Compute intersection: */
                if ((aux->t.zMin < planes[p]) && (aux->t.zMax > planes[p])) {
                    LineSegment seg = R3_Mesh_Triangle_slice(aux, planes[p]);
                    segs[p].push_back(seg);
//                    intersections++;
                }
            }
            aux = next;
        }
    }
    free(L);
    clock_t slice_end = clock();
    double slicing_time = double(slice_end - slice_begin) / CLOCKS_PER_SEC;
    std::cout << "Incremental " << slicing_time << std::endl;

    /*End-Slicing*/

//    if (chaining) {
//        /*Contour construction:*/
//        for (size_t p = 0; p < k; p++) {
//            if (!segs[p].empty()) {
//                ContourConstruction(segs[p], polygons, p);
//                if (orienting) {
//                    ray_casting(polygons[p]);
//                }
//                segs[p].clear();
//            }
//        }
//        /*End construction.*/
//    } else {
//        export_svg_no_chaining("segments.svg", segs, k, mesh->meshAABBSize());
//    }
    std::vector<Layer *> result;

    return result;
}


std::vector<Layer *> PolygonSlicer::sliceModel(TriangleMesh& mesh, double thickness, double epsilon) {
    std::vector<float> planes = computePlanes(mesh, thickness, epsilon);

    Timer timer;
    std::vector<Layer *> result;

    result = TrivialSlicing(mesh, planes);
//    std::cout << "Trivial Time: " << timer.elapsed() << " seconds" << std::endl;

//    timer.reset();
//    result = IncrementalSlicing(mesh, planes, thickness);
//    std::cout << "Incremental Time: " << timer.elapsed() << " seconds" << std::endl;
    return result;
}
