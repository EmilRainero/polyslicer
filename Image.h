//
// Created by erainero on 11/20/18.
//

#ifndef RASTERIZE_IMAGE_H
#define RASTERIZE_IMAGE_H

#include "rasterize/Color.h"

class Image {
public:
    Image(int width, int height);

    void setPixels(int y, int x1, int x2, const Color& color);

    void print();

    int width;
    int height;
    unsigned char* buffer;
};


#endif //RASTERIZE_IMAGE_H
