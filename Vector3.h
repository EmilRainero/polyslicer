
#ifndef POLYSLICER_VECTOR3_H
#define POLYSLICER_VECTOR3_H

#include <iostream>
#include <sstream>
#include <cmath>

class Vector3 {

public:

    Vector3(float _x = 0, float _y = 0, float _z = 0) : x(_x), y(_y), z(_z) {}

    float distTo(const Vector3 &pt) const {
        return sqrt(pow(fabs(x - pt.x), 2.0) + pow(fabs(y - pt.y), 2.0) + pow(fabs(z - pt.z), 2.0));
    }

    Vector3 &operator-=(const Vector3 &pt) {
        x = (x - pt.x);
        y = (y - pt.y);
        z = (z - pt.z);
        return *this;
    }

    bool operator==(const Vector3 &pt) const {
        return distTo(pt) < 0.005;
    }

    std::string getLabel() const {
        std::stringstream ss;
        ss << x << "|" << y << "|" << z;
        return ss.str();
    }

    friend std::ostream &operator<<(std::ostream &os, const Vector3 &v) {
        os << "x: " << v.x << "; y: " << v.y << "; z: " << v.z;
        return os;
    }

public:

    float x;
    float y;
    float z;
};

#endif //POLYSLICER_VECTOR3_H
