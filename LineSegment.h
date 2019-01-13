#ifndef POLYSLICER_LINESEGMENT_H
#define POLYSLICER_LINESEGMENT_H

#include "Vector3.h"

class LineSegment {

public:

    LineSegment (Vector3 p0=Vector3(), Vector3 p1=Vector3(), int i=0) {
        v[0] = p0;
        v[1] = p1;
        index = i;
        vertical = false;
        if ((v[1].x - v[0].x) != 0) {
            a = (v[1].y - v[0].y)/(v[1].x - v[0].x);
            b = (v[0].y - (a * v[0].x));
        } else {
            vertical = true;
        }
    }

    bool operator==(const LineSegment &ls) const {
        return ((v[0] == ls.v[0]) && (v[1] == ls.v[1]));
    }

    friend std::ostream& operator<<(std::ostream& os, const LineSegment& ls) {
        os << "V0: (" << ls.v[0] << "); V1: (" << ls.v[1] << ")";
        return os;
    }

public:

    Vector3 v[2];
    double a;
    double b;
    bool vertical;
    int index;
};

#endif //POLYSLICER_LINESEGMENT_H
