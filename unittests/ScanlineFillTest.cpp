#include "gtest/gtest.h"
#include "ScanlineFill.h"
#include "Point.h"
#include "Polygon.h"
#include "Image.h"


class ScanlineFillTest : public ::testing::Test {
protected:
    ScanlineFillTest() : backgroundColor(0), foregroundColor(255) {

    }

    ~ScanlineFillTest() override = default;

    Color backgroundColor;
    Color foregroundColor;

    void createRectangle(Polygon2i& polygon, Point2i lowerLeft, Point2i upperRight) {
        polygon.pt.push_back(Point2i(lowerLeft.x, lowerLeft.y));
        polygon.pt.push_back(Point2i(lowerLeft.x, upperRight.y));
        polygon.pt.push_back(Point2i(upperRight.x, upperRight.y));
        polygon.pt.push_back(Point2i(upperRight.x, lowerLeft.y));
        polygon.pt.push_back(Point2i(lowerLeft.x, lowerLeft.y));
    }
};




TEST_F(ScanlineFillTest, Test1) {
    int imageSize{30};
    int padding{10};
    Image image(imageSize, imageSize);
    Point2i lowerLeft(padding, padding);
    Point2i upperRight(imageSize-padding, imageSize-padding);

    Polygon2i polygon;
    createRectangle(polygon, lowerLeft, upperRight);

    ScanlineFill::scanFill(polygon, image, foregroundColor);
//    image.print();

    for (int x = 0; x < image.width; x++) {
        for (int y = 0; y < image.height; y++) {
            bool pixelShouldBeSet = (lowerLeft.x <= x && x < upperRight.x ) &&
                                    (lowerLeft.y <= y && y < upperRight.y);
            ASSERT_EQ(pixelShouldBeSet ? 255 : 0, image.getPixel(x, y));
        }
    }
}

TEST_F(ScanlineFillTest, Test2) {
    int imageSize{10000};
    int padding{10};
    Image image(imageSize, imageSize);
    Point2i lowerLeft(padding, padding);
    Point2i upperRight(imageSize-padding, imageSize-padding);

    Polygon2i polygon;
    polygon.pt.push_back(Point2i(lowerLeft.x, lowerLeft.y));
    for (int x = lowerLeft.x; x < upperRight.x; x+= 20) {
    polygon.pt.push_back(Point2i(x, upperRight.y));
    polygon.pt.push_back(Point2i(x+10, upperRight.y));
    polygon.pt.push_back(Point2i(x+10,lowerLeft.y+10));
    polygon.pt.push_back(Point2i(x+20,lowerLeft.y+10));
    }
    polygon.pt.push_back(Point2i(upperRight.x,lowerLeft.y));
    polygon.pt.push_back(Point2i(lowerLeft.x, lowerLeft.y));

    ScanlineFill::scanFill(polygon, image, foregroundColor);
//        image.print();
}
