//
// Created by erainero on 1/12/19.
//

#ifndef POLYSLICER_BOUNDINGBOX_H
#define POLYSLICER_BOUNDINGBOX_H

class BoundingBox {
public:
    BoundingBox() {
        xMax = std::numeric_limits<double>::min();
        xMin = std::numeric_limits<double>::max();
        yMax = std::numeric_limits<double>::min();
        yMin = std::numeric_limits<double>::max();
    }

    double xMin;
    double xMax;
    double yMin;
    double yMax;
    bool firstPoint;
};

#endif //POLYSLICER_BOUNDINGBOX_H
