
#ifndef POLYSLICER_VECTOR3_H
#define POLYSLICER_VECTOR3_H

#include <iostream>
#include <sstream>
#include <cmath>

class Vector3 {

public:

    Vector3 (float _x=0, float _y=0, float _z=0) : x(_x), y(_y), z(_z) {}

    float distTo (const Vector3 &pt) const {
        return sqrt ( pow(fabs(x-pt.x), 2.0) + pow(fabs(y-pt.y), 2.0) + pow(fabs(z-pt.z), 2.0) );
    }

//    array<float,3> getCoords() {
//        array<float,3> c = {{x, y, z}};
//        return c;
//    }

    float dotproduct (const Vector3 &v) const {
        return (x*v.x + y*v.y + z*v.z);
    }

//    void transform (const glm::mat4 &mat) {
//        glm::vec4 v = glm::vec4(x, y, z, 1.0);
//        glm::vec4 vt = mat*v;
//        x = (vt.x); y = (vt.y); z = (vt.z);
//    }

    Vector3& operator-=(const Vector3 &pt) {
        x = (x-pt.x);
        y = (y-pt.y);
        z = (z-pt.z);
        return *this;
    }

    Vector3 operator-(const Vector3 &pt) {
        return Vector3 ((x-pt.x), (y-pt.y), (z-pt.z));
    }

    Vector3 operator+(const Vector3 &pt) {
        return Vector3 ((x+pt.x), (y+pt.y), (z+pt.z));
    }

    Vector3 operator/(float a) {
        return Vector3 ((x/a), (y/a), (z/a));
    }

    Vector3 operator*(float a) {
        return Vector3 ((x*a), (y*a), (z*a));
    }

    bool operator<(const Vector3 &pt) const {
        return z < pt.z;
    }

    bool operator>(const Vector3 &pt) const {
        return z > pt.z;
    }

    bool operator==(const Vector3 &pt) const {
        return distTo(pt) < 0.005;
    }

    bool operator!=(const Vector3 &pt) const {
        return distTo(pt) > 0.005;
    }

    float normalize() const {
        return sqrt(x*x+y*y+z*z);
    }

    std::string getLabel() const {
        std::stringstream ss;
        ss << x << "|" << y << "|" << z;
        return ss.str();
    }

    friend std::ostream& operator<<(std::ostream& os, const Vector3& v) {
        os << "x: " << v.x << "; y: " << v.y << "; z: " << v.z;
        return os;
    }

public:

    float x;
    float y;
    float z;
};

#endif //POLYSLICER_VECTOR3_H
