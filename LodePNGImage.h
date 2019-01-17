#ifndef POLYSLICER_LODEPNGIMAGE_H
#define POLYSLICER_LODEPNGIMAGE_H

#include <cmath>
#include <vector>
#include <string>

class LodePNGImage {
private:
    double pixelDifference = 1.0;

public:
    std::vector<unsigned char> data;
    unsigned width, height;
    unsigned size;
    std::string filename;

    LodePNGImage(unsigned _width, unsigned _height, int value = 255);

    LodePNGImage(const std::string _filename);

    ~LodePNGImage();

    void write(const char *filename);

    void setPixel(int x, int y, int r, int g, int b);

    double computeDifference(double r1, double g1, double b1, double r2, double g2, double b2);

    bool deepCompare(const LodePNGImage& image2, const double differencePercentage, std::string& error);

    LodePNGImage difference(const LodePNGImage& image2);

    int arrayIndex(const int x, const int y);

    void erode();

    void dilate();

    int countPixelsOfColor(int r, int g, int b);
};

#endif //POLYSLICER_LODEPNGIMAGE_H
