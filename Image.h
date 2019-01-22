//
// Created by erainero on 11/20/18.
//

#ifndef RASTERIZE_IMAGE_H
#define RASTERIZE_IMAGE_H

#include "Color.h"

class Image {
public:
    Image(int width, int height);

    ~Image();

    void setPixels(int y, int x1, int x2, const Color& color);
    int getPixel(int x, int y) {
        return buffer[y*width + x];
    }

    void print();

    int width;
    int height;
    unsigned char* buffer;
};


#endif //RASTERIZE_IMAGE_H
