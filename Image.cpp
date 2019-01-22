//
// Created by erainero on 11/20/18.
//

#include <iostream>
#include <string.h>
#include <stdlib.h>

#include "Image.h"

Image::Image(int width, int height) : width(width), height(height) {
    buffer = (unsigned char*) calloc(height * width, sizeof(unsigned char));
}

Image::~Image() {
    free(buffer);
}

void Image::setPixels(int y, int x1, int x2, const Color& color) {
//    std::cout << y << "  " << x1 << "-" << x2 << std::endl;
    if (    y < 0 || y >= height ||
                     x1  < 0 || x1 >= width ||
                                x2 < 0 || x2 >= width
            ) {
        std::cerr << "setPixels error " << y << "  " << x1 << "-" << x2 << std::endl;
        return;
    }
    if (x1 > x2) {
        int swap = x1;
        x1 = x2;
        x2 = swap;
    }
    memset(buffer+y*width + x1, color.color, x2-x1);
}

void Image::print() {
    for (int y = 0; y < height; y++) {
        for (int x = 0; x< width; x++) {
            std::cout << (buffer[y*width + x] == 0 ? ' ' : 'X');
        }
        std::cout << std::endl;
    }
}
