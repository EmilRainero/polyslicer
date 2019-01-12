//
// Created by erainero on 1/12/19.
//

#ifndef POLYSLICER_CONTOUR_H
#define POLYSLICER_CONTOUR_H

class Contour {
public:
    bool external;
    bool clockwise;
    std::vector<Vector3> points;
};

#endif //POLYSLICER_CONTOUR_H
