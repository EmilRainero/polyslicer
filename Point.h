//
// Created by erainero on 11/20/18.
//

#ifndef RASTERIZE_POINT2I_H
#define RASTERIZE_POINT2I_H

#include <iostream>
#include <fstream>


template <class T>
class Point2 {
public:
    T x, y;    // elements

    Point2() {}
    Point2(T x, T y) : x(x), y(y) {}

    void writeBinary(std::ofstream& fileStream) {
        fileStream.write(reinterpret_cast<char*>(&x), sizeof(x));
        fileStream.write(reinterpret_cast<char*>(&y), sizeof(y));
    }

    static Point2<T> readBinary(std::ifstream& fileStream) {
        Point2<T> point;

        fileStream.read(reinterpret_cast<char*>(&point.x), sizeof(point.x));
        fileStream.read(reinterpret_cast<char*>(&point.y), sizeof(point.y));
        return point;
    }

    void print() {
        std::cout << "(" << x << ", " << y << ") ";
    }

};

typedef Point2<int> Point2i;
typedef Point2<float> Point2f;
typedef Point2<double> Point2d;

template <class T>
class BBox2 {
public:
    T min;
    T max;

    BBox2() {}

    BBox2(T min, T max) : min(min), max(max) {}
    void writeBinary(std::ofstream& fileStream) {
        std::cout << "Point3 writeBinary" << std::endl;
//        fileStream.write(reinterpret_cast<char*>(&numberPolygons), sizeof(numberPolygons));
    }

};

typedef BBox2<int> BBox2i;

template <class T>
class Point3 {
public:
    T x, y, z;    // elements

public:
    Point3() {}
    Point3(T x, T y, T z) : x(x), y(y), z(z) {}

    void print() {
        std::cout << "(" << x << ", " << y << ", " << z << ") ";
    }
};

typedef Point3<int> Point3i;
typedef Point3<float> Point3f;
typedef Point3<double> Point3d;

#endif //RASTERIZE_POINT2I_H
