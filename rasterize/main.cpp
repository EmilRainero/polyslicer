// 2D Scan Line Fill Using an Edge Table and Active Edge List
// One of the standard scan line fill algorithms uses an Edge Table and Active Edge List. We build these data structures using the C++ STL containers of vector and list. The code follows.

#include <iostream>
#include <cmath>
#include <iomanip>

#include "../Point.h"
#include "../Polygon.h"
#include "../Node.h"
#include "../EdgeTable.h"
#include "../Image.h"

#include "JsonUtils.h"
#include "Timer.h"



void CreatePolygonFromJSON(Polygon2i &P, Json::Value& points) {
    Point2i p;
    Point2i lastp;

    P.pt.clear();
    for (int j = 0; j < points.size(); j++) {
        Json::Value point1 = points[j];
        p.x = round(point1[0].asDouble());
        p.y = round(point1[1].asDouble());
        if (j > 0 && !(lastp.x == p.x && lastp.y == p.y)) {
            P.pt.push_back(p);
        }
        lastp = p;
    }

    // remove horizontal or vertical redundant points
    for (int j = P.pt.size() - 2; j > 0; j--) {
        if (
                (P.pt[j - 1].x == P.pt[j].x && P.pt[j].x == P.pt[j + 1].x) ||
                (P.pt[j - 1].y == P.pt[j].y && P.pt[j].y == P.pt[j + 1].y)
                )
            P.pt.erase(P.pt.begin() + j);
    }
}

void printPolygonArray(PolygonArray2i& polygonArray) {
    std::cout << "Number polygons: " << polygonArray.polygons.size() << std::endl;
    for (auto polygon: polygonArray.polygons) {
        std::cout << "\t" << polygon.pt.size() << " - ";
        for (auto point: polygon.pt) {
            std::cout << point.x << "," << point.y << " ";
        }
        std::cout << std::endl;
    }
}

void testBinaryReadyWrite(std::string filename) {
    Timer timer;
    timer.reset();

    std::cout << std::fixed;

    PolygonArray2i polygonArray = PolygonArray2i::readJSON(filename);

    std::cout << std::setprecision(10) << "Read JSON Time:    " << timer.elapsed() << std::endl;

    timer.reset();
    polygonArray.writeBinary("polygons.bin");
    std::cout << std::setprecision(10) << "Write Binary Time: " << timer.elapsed() << std::endl;

    timer.reset();
    PolygonArray2i tempPolygonArray = PolygonArray2i::readBinary("polygons.bin");
    std::cout << std::setprecision(10) << "Rad Binary Time:   " << timer.elapsed() << std::endl;
//    printPolygonArray(tempPolygonArray);

}

void runBinary(std::string filename, int width, int height) {
    Image image(width, height);

    Color backgroundColor(0);
    Color foregroundColor(255);

    Timer timer;
    timer.reset();

    PolygonArray2i polygonArray = PolygonArray2i::readBinary(filename);

    for (int i = 0; i< polygonArray.polygons.size(); i++) {
        Polygon2i polygon = polygonArray.polygons[i];
//        std::cout << "   " << polygon.pt.size() << std::endl;

        EdgeTable::scanFill(polygon, image, polygon.interior ? backgroundColor : foregroundColor);
    }

    std::cout << "Time: " << timer.elapsed() << std::endl;
//    image.print();

//    printPolygonArray(polygonArray);
//
//    polygonArray.writeBinary("polygons.bin");
//    PolygonArray2i tempPolygonArray = PolygonArray2i::readBinary("polygons.bin");
//
//    printPolygonArray(tempPolygonArray);

}

void run(std::string filename, int width, int height) {
    Image image(width, height);

    Color backgroundColor(0);
    Color foregroundColor(255);

    Timer timer;
    timer.reset();

    PolygonArray2i emil;

    PolygonArray2i polygonArray = PolygonArray2i::readJSON(filename);

    for (int i = 0; i< polygonArray.polygons.size(); i++) {
        Polygon2i polygon = polygonArray.polygons[i];

        EdgeTable::scanFill(polygon, image, polygon.interior ? backgroundColor : foregroundColor);
    }

    std::cout << "Time: " << timer.elapsed() << std::endl;
    image.print();

}

void convertJsonToBinary(std::string jsonFilename, std::string binaryFilename) {
    PolygonArray2f polygonArray = PolygonArray2f::readJSON(jsonFilename);

    polygonArray.print();
    auto bbox = polygonArray.computeBoundingBox();
    std::cout <<
              bbox.min.x << ", " << bbox.min.y << "  " <<
              bbox.max.x << ", " << bbox.max.y << std::endl;
}

void boundingBox(std::string filename, int width, int height) {
    Image image(width, height);

    PolygonArray2i polygonArray = PolygonArray2i::readBinary(filename);
    auto bbox = polygonArray.computeBoundingBox();
    std::cout <<
              bbox.min.x << ", " << bbox.min.y << "  " <<
              bbox.max.x << ", " << bbox.max.y << std::endl;
}

int main(int argc, char **argv) {
//    convertJsonToBinary("../testbox.json", "testbox.bin");
//    convertJsonToBinary("../00076.json", "00076.bin");

//    boundingBox("00076.bin", 1024, 1024);

//    return 0;

    std::cout << std::fixed;

//    runBinary("testbox.bin", 11000, 11000);
//    testBinaryReadyWrite("../00076.json");

    run("../00076.json", 1024, 1024);
//    run("../engine/00009.json", 1024, 1024);
//    run("../testbox.json", 11000, 11000);

    return 0;
}


