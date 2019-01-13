
#ifndef POLYSLICER_TRIANGLE_H
#define POLYSLICER_TRIANGLE_H

#include "Vector3.h"

class Triangle    {

public:

    Triangle(Vector3 n, Vector3 v0, Vector3 v1, Vector3 v2) : normal(n) {
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

    Triangle& operator-=(const Vector3 &pt) {
        v[0] -= pt;
        v[1] -= pt;
        v[2] -= pt;
        return *this;
    }

    bool operator<(const Triangle &t) {
        return zMin < t.zMin;
    }

    friend std::ostream& operator<<(std::ostream& os, const Triangle& t) {
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

    Vector3 v[3];
    Vector3 normal;
    float zMin;
    float zMax;
};

#endif
