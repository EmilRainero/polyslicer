//
// Created by erainero on 11/20/18.
//

#ifndef RASTERIZE_POLYGON2I_H
#define RASTERIZE_POLYGON2I_H

#include <vector>
#include <string>
#include <fstream>
#include <cmath>

#include "Point.h"
#include "JsonUtils.h"


template <class P>
class Polygon {
public:
    std::vector<P> pt;

    bool interior;

    void writeBinary(std::ofstream& fileStream);

    static Polygon<P> readBinary(std::ifstream& fileStream);

    void print();
};

template <class P>
void Polygon<P>::print() {
    std::cout << "Polygon " << pt.size() << " - ";
    for (P point: pt) {
        point.print();
    }
    std::cout << std::endl;
}

template <class P>
void Polygon<P>::writeBinary(std::ofstream& fileStream) {
//    std::cout << "  polygon. writeBinary" << std::endl;

    int numberPoints = pt.size();
//    std::cout << "write points " << numberPoints << std::endl;
    fileStream.write(reinterpret_cast<char*>(&numberPoints), sizeof(numberPoints));
    for (auto point: pt) {
        point.writeBinary(fileStream);
    }
}

template <class P>
Polygon<P> Polygon<P>::readBinary(std::ifstream& fileStream) {
    Polygon<P> polygon;

    int numberPoints;
    fileStream.read(reinterpret_cast<char*>(&numberPoints), sizeof(numberPoints));
//    std::cout << "read points " << numberPoints << std::endl;
    for (int i = 0; i< numberPoints; i++) {
        P pt = P::readBinary(fileStream);
        polygon.pt.push_back(pt);
    }
    return polygon;
}

typedef Polygon<Point2i> Polygon2i;
typedef Polygon<Point2f> Polygon2f;


template <class P, class PT>
class PolygonArray {
public:
    std::vector<P> polygons;

    void print();

    static PolygonArray<P, PT> readJSON(std::string filename);

    P CreatePolygonFromJSON(Json::Value& points);

    void writeBinary(std::string filename);

    static PolygonArray<P, PT> readBinary(std::string filename);

    BBox2<PT> computeBoundingBox();
};

template <class P, class PT>
BBox2<PT> PolygonArray<P, PT>::computeBoundingBox() {
    BBox2<PT> bbox;
    bbox.min = polygons[0].pt[0];
    bbox.max = polygons[0].pt[0];
    for (auto polygon: polygons) {
        for (auto point: polygon.pt) {
            bbox.min.x = std::min(bbox.min.x, point.x);
            bbox.min.y = std::min(bbox.min.y, point.y);
            bbox.max.x = std::max(bbox.max.x, point.x);
            bbox.max.y = std::max(bbox.max.y, point.y);
        }
    }
    return bbox;
}

template <class P, class PT>
void PolygonArray<P, PT>::print() {
    for (auto polygon: polygons) {
        polygon.print();
    }
}

template <class P, class PT>
PolygonArray<P, PT> PolygonArray<P, PT>::readJSON(std::string filename) {
    Json::Value polygonsJSON = JsonUtils::readFromFile(filename);

    PolygonArray<P, PT> polygonArray;

    for (int i = 0; i < polygonsJSON.size(); i++) {
        Json::Value polygonJSON = polygonsJSON[i];
        Json::Value pointsJSON = polygonJSON["Points"];
        bool interior = polygonJSON["Interior"].asBool();

        Polygon<PT> polygon;
        PT p;
        PT lastp;

        for (int j = 0; j < pointsJSON.size(); j++) {
            Json::Value point1 = pointsJSON[j];
            p.x = (point1[0].asDouble());
            p.y = (point1[1].asDouble());
            if (j > 0 && !(lastp.x == p.x && lastp.y == p.y)) {
                polygon.pt.push_back(p);
            }
            lastp = p;
        }
        // remove horizontal or vertical redundant points
        for (int j = polygon.pt.size() - 2; j > 0; j--) {
            if (
                    (polygon.pt[j - 1].x == polygon.pt[j].x && polygon.pt[j].x == polygon.pt[j + 1].x) ||
                    (polygon.pt[j - 1].y == polygon.pt[j].y && polygon.pt[j].y == polygon.pt[j + 1].y)
                    )
                polygon.pt.erase(polygon.pt.begin() + j);
        }

        polygon.interior = interior;
        polygonArray.polygons.push_back(polygon);
    }

    return polygonArray;
}

template <class P, class PT>
void PolygonArray<P, PT>::writeBinary(std::string filename) {
//    std::cout << "polygonarray writeBinary" << std::endl;
    std::ofstream fileStream(filename, std::ios::out | std::ios::binary);

    int numberPolygons = polygons.size();
    fileStream.write(reinterpret_cast<char*>(&numberPolygons), sizeof(numberPolygons));
    for (auto polygon: polygons) {
        polygon.writeBinary(fileStream);
    }
    fileStream.close();
}

template <class P, class PT>
PolygonArray<P, PT> PolygonArray<P, PT>::readBinary(std::string filename) {
    PolygonArray<P, PT> polygonArray;

    std::ifstream fileStream(filename, std::ios::in | std::ios::binary);

    int numberPolygons;
    fileStream.read(reinterpret_cast<char*>(&numberPolygons), sizeof(numberPolygons));

    for (int i = 0; i< numberPolygons; i++) {
        Polygon<PT> polygon = Polygon<PT>::readBinary(fileStream);
        polygonArray.polygons.push_back(polygon);
    }
    fileStream.close();

    return polygonArray;
}



typedef PolygonArray<Polygon2i, Point2i> PolygonArray2i;
typedef PolygonArray<Polygon2f, Point2f> PolygonArray2f;


#endif //RASTERIZE_POLYGON2I_H
