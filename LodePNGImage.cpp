#include <algorithm>
#include <stdexcept>
#include "LodePNGImage.h"
#include "LodePNG.h"

LodePNGImage::LodePNGImage(unsigned _width, unsigned _height, int value) {
    width = _width;
    height = _height;
    size = width * height;
    data.resize(width * height * 4, value); // white

    // alpha
    for (unsigned i = 0; i < width * height; i++) {
        data[4 * i + 3] = 255;
    }
}

LodePNGImage::LodePNGImage(const std::string _filename) {
    unsigned error;

    filename = _filename;
    error = lodepng::decode(data, width, height, _filename);
    if (error) {
        throw std::invalid_argument("Error " + std::to_string(error));
    }
}

LodePNGImage::~LodePNGImage() = default;

double LodePNGImage::computeDifference(double r1, double g1, double b1, double r2, double g2, double b2) {

    double rdiff = fabs(r1 / 255.0 - r2 / 255.0);
    double gdiff = fabs(g1 / 255.0 - g2 / 255.0);
    double bdiff = fabs(b1 / 255.0 - b2 / 255.0);

    double diff = rdiff + gdiff + bdiff;
    diff = fmin(1.0, diff * 2);

    return diff;
}

void LodePNGImage::write(const char *filename) {
    lodepng::encode(filename, data, width, height);
}

void LodePNGImage::setPixel(int x, int y, int r, int g, int b) {
    int index = arrayIndex(y, x);
    data[index] = r;
    data[index+1] = g;
    data[index+2] = b;
}

bool LodePNGImage::deepCompare(const LodePNGImage& image2, const double differencePercentage, std::string& error) {
    LodePNGImage result = difference(image2);

    result.erode();

//    lodepng::encode("/tmp/out.png", result.data, result.width, result.height);

    int count = count = result.countPixelsOfColor(0, 0, 0);
    double actualDifferencePercentage = ((double) count) / (double) result.size;

    bool success = actualDifferencePercentage <= differencePercentage;

    if (!success) {
        error = "files differ " + filename + " " + image2.filename + " percentage different: " +
                std::to_string(actualDifferencePercentage * 100.0);
    }

    return success;
}

LodePNGImage LodePNGImage::difference(const LodePNGImage& image2) {
    unsigned size = width * height;

    LodePNGImage result(width, height);

    for (unsigned i = 0; i < size; i++) {
        int r1, g1, b1, a1;
        int r2, g2, b2, a2;

        r1 = data[4 * i];
        g1 = data[4 * i + 1];
        b1 = data[4 * i + 2];
        a1 = data[4 * i + 3];
        r2 = image2.data[4 * i];
        g2 = image2.data[4 * i + 1];
        b2 = image2.data[4 * i + 2];
        a2 = image2.data[4 * i + 3];

        double difference = computeDifference((double) r1, (double) g1, (double) b1, (double) r2, (double) g2,
                                              (double) b2);
        if (difference >= pixelDifference) {
            result.data[4 * i] = 0;
            result.data[4 * i + 1] = 0;
            result.data[4 * i + 2] = 0;
        }
    }

    return result;
}

int LodePNGImage::arrayIndex(const int x, const int y) {
    return (y * width + x) * 4;
}

void LodePNGImage::erode() {

    unsigned size = width * height;
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            if (data[arrayIndex(i, j)] == 0) {
                if (i > 0 && data[arrayIndex(i - 1, j)] == 255)
                    data[arrayIndex(i - 1, j)] = 1;
                if (j > 0 && data[arrayIndex(i, j - 1)] == 255)
                    data[arrayIndex(i, j - 1)] = 1;
                if (i + 1 < width && data[arrayIndex(i + 1, j)] == 255)
                    data[arrayIndex(i, j)] = 1;
                if (j + 1 < height && data[arrayIndex(i, j + 1)] == 255)
                    data[arrayIndex(i, j + 1)] = 1;
            }
        }
    }
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (data[arrayIndex(i, j)] == 1) {
                data[arrayIndex(i, j)] = 255;
                data[arrayIndex(i, j) + 1] = 255;
                data[arrayIndex(i, j) + 2] = 255;
            }
        }
    }
}

void LodePNGImage::dilate() {
    unsigned size = width * height;
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            if (data[arrayIndex(i, j)] == 0) {
                if (i > 0 && data[arrayIndex(i - 1, j)] == 255)
                    data[arrayIndex(i - 1, j)] = 1;
                if (j > 0 && data[arrayIndex(i, j - 1)] == 255)
                    data[arrayIndex(i, j - 1)] = 1;
                if (i + 1 < width && data[arrayIndex(i + 1, j)] == 255)
                    data[arrayIndex(i, j)] = 1;
                if (j + 1 < height && data[arrayIndex(i, j + 1)] == 255)
                    data[arrayIndex(i, j + 1)] = 1;
            }
        }
    }
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (data[arrayIndex(i, j)] == 1) {
                data[arrayIndex(i, j)] = 0;
                data[arrayIndex(i, j) + 1] = 0;
                data[arrayIndex(i, j) + 2] = 0;
            }
        }
    }
}

int LodePNGImage::countPixelsOfColor(int r, int g, int b) {
    int count = 0;
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            int index = arrayIndex(i, j);
            if (data[index] == r && data[index + 1] == g && data[index + 2] == b) {
                count++;
            }
        }
    }
    return count;
}
