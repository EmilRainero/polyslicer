#include "gtest/gtest.h"
#include "EdgeTable.h"
#include "Point.h"
#include "Polygon.h"
#include "Image.h"


class EdgeListTest : public ::testing::Test {
    protected:
    EdgeListTest() : backgroundColor(0), foregroundColor(255) {

    }

    ~EdgeListTest() override = default;

    Color backgroundColor;
    Color foregroundColor;
};

TEST_F(EdgeListTest, Test1) {
    Image image(30, 30);
    Point2i lowerLeft(10, 10);
    Point2i upperRight(20, 20);

    Polygon2i polygon;
    polygon.pt.push_back(Point2i(lowerLeft.x, lowerLeft.y));
    polygon.pt.push_back(Point2i(lowerLeft.x, upperRight.y));
    polygon.pt.push_back(Point2i(upperRight.x, upperRight.y));
    polygon.pt.push_back(Point2i(upperRight.x, lowerLeft.y));
    polygon.pt.push_back(Point2i(lowerLeft.x, lowerLeft.y));

    EdgeTable::scanFill(polygon, image, foregroundColor);

//    image.print();

    for (int x = 0; x < image.width; x++) {
        for (int y = 0; y < image.height; y++) {
            bool pixelShouldBeSet = (lowerLeft.x <= x && x < upperRight.x ) &&
                    (lowerLeft.y <= y && y < upperRight.y);
            ASSERT_EQ(pixelShouldBeSet ? 255 : 0, image.getPixel(x, y));
        }
    }
    ASSERT_TRUE(true);


}
