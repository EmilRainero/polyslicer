//
// Created by erainero on 1/12/19.
//

#ifndef POLYSLICER_LAYER_H
#define POLYSLICER_LAYER_H

#include "Contour.h"

class Layer {
public:

    Layer(float z) : z(z) {}

    float z;
    std::vector<Contour> contours;
};

#endif //POLYSLICER_LAYER_H
