#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>


#include "PolygonSlicer.h"
#include "Timer.h"
#include "Image.h"
#include "EdgeTable.h"
#include "LodePNG.h"
#include "LodePNGImage.h"



std::vector<Vector3> minimizeContour(Contour& contour) {

    std::vector<Vector3>& points = contour.points;

    int size = points.size() - 1; // ignore the extra point at end

    std::vector<Vector3>newPoints;
    for (int i = 0; i< size; i++) {
        if (i > 0 && points[i].x == points[i-1].x && points[i].y == points[i-1].y) {
            continue;
        }
        int previousIndex = (i+size-1) % size;
        int nextIndex = (i+1) % size;

//        std::cout << "trying "
//                << points[previousIndex].x << "," << points[previousIndex].y << "   "
//                << points[i].x << "," << points[i].y << "   "
//                << points[nextIndex].x << "," << points[nextIndex].y << " "
//                << points[i].isColinear(points[previousIndex], points[nextIndex])
//                     << std::endl;
        if (!points[i].isColinear(points[previousIndex], points[nextIndex])) {
            newPoints.push_back(points[i]);
//            std::cout << points[i].x << "," << points[i].y << std::endl;
        }
    }
    newPoints.push_back(newPoints[0]);
    std::cout << "Before " << size+1 << "  After " << newPoints.size() << std::endl;

    return newPoints;
}

void print(Contour& contour) {
    std::cout << "\tContour - external: " << contour.external << " - clockwise: " << contour.clockwise
              << std::endl;

//    for (auto point: contour.points) {
//        std::cout << "\t\t" << point << std::endl;
//    }

    auto data = minimizeContour(contour);
//    for (auto point: data) {
//        std::cout << "\t\t" << point << std::endl;
//    }
}

void print(Layer* layer) {
    for (auto contour: layer->contours) {
        print(contour);
    }
}

void print(std::vector<Layer *> layers) {
    for (auto layer: layers) {
        std::cout << "Layer: " << layer->z << std::endl;

        print(layer);
    }
}

void rescaleLayer(Layer *layer, float tx, float ty, float scale) {
    for (Contour &contour: layer->contours) {
        for (Vector3 &point: contour.points) {
            point.x = (point.x - tx) * scale;
            point.y = (point.y - ty) * scale;
        }
    }
}

void simplifyPoints(std::vector<P>& points) {

}

Image rasterizeLayer(Layer *layer, int xSize, int ySize) {
    Image image(xSize, ySize);

    Color backgroundColor(0);
    Color foregroundColor(255);

//    std::cout << "Rasterize " << layer->z << std::endl;

    for (int i = 0; i < 2; i++) {

        for (Contour &contour: layer->contours) {

            if (
                    (i == 0 && !contour.external) ||
                    (i == 1 && contour.external)) {
                continue;
            }
            Polygon2i polygon;

            for (Vector3 &point: contour.points) {
                Point2i pt(point.x, point.y);
                polygon.pt.push_back(pt);
            }
            simplifyPoints(polygon.pt);

            polygon.interior = !contour.external;

//            typedef struct {
//                int x, y;
//            } Point;
//            int data[][2] = {
//                    {0,1},
//                    {1,11},
//                    {11,11},
//                    {11,1},
//                    {0,1},
//            };
//            int n = sizeof(data) / sizeof(int *);
//            polygon.pt.clear();
//            for (int i = 0; i < n; i++) {
//                Point2i p(data[i][0], data[i][1]);
//                polygon.pt.push_back((p));
//            }

            /*
             *   when 1
             *   contents: -1 0 0
             *   contents: 11 1 0
             *   contents: 11 11 0
             *
             *   when 0
             *   contents: 11 0 0.1
             *   contents: -1 0 0
             *   contents: 11 11 0
             */

//            std::cout << "fill " << polygon.interior << std::endl;
            EdgeTable::scanFill(polygon, image, polygon.interior ? backgroundColor : foregroundColor);
        }
    }

    return image;
}

std::string ZeroPadNumber(int number, int digits) {
    std::ostringstream ss;
    ss << std::setw(digits) << std::setfill('0') << number;
    return ss.str();
}

void rasterizeLayers(TriangleMesh &mesh, std::vector<Layer *> layers) {
    int arraySize{1024};

    auto minVertex = mesh.getBottomLeftVertex();
    auto maxVertex = mesh.getTopRightVertex();

    float xSize = maxVertex.x - minVertex.x;
    float ySize = maxVertex.y - minVertex.y;
    float maxScale = std::max(xSize, ySize);

//    std::cout << xSize << " " << ySize << " " << maxScale << std::endl;

    float scale = (arraySize - 1) / maxScale;

    int layerNumber = 0;

    double totalRasterizeTime{0.0};
    double totalWriteTime{0.0};
    for (auto layer: layers) {
        Timer timer;

        std::string filename = "layer-" + ZeroPadNumber(layerNumber, 5) + ".png";
        std::cout << filename;

        rescaleLayer(layer, minVertex.x, minVertex.y, scale);
        auto image = rasterizeLayer(layer, arraySize, arraySize);
        double elapsed;
        elapsed = timer.elapsed();
        totalRasterizeTime += elapsed;
        std::cout << " " << elapsed;
        if (true) {
            timer.reset();
            LodePNGImage pngImage(arraySize, arraySize, 0);
            for (int y = 0; y < arraySize; y++) {
                for (int x = 0; x < arraySize; x++) {
                    if (image.buffer[y * arraySize + x] > 0) {
                        pngImage.setPixel(x, y, 255, 255, 255);
                    }
                }
            }

            pngImage.write(filename.c_str());
            elapsed = timer.elapsed();
            totalWriteTime += elapsed;

            std::cout << "  Write " << elapsed;
        }

        std::cout << std::endl;

//        image.print();
        layerNumber++;
    }

    std::cout << "Rasterize time: " << totalRasterizeTime << "  Write time: " << totalWriteTime << std::endl;
}

int main() {
    double epsilon{0.0001};
    double layerHeight{1};
    Timer timer;

//    TriangleMesh mesh("../parts/cube_10x10x10.stl", epsilon);
//    TriangleMesh mesh("../parts/bblocky.stl", epsilon);
    TriangleMesh mesh("../parts/pisa.stl", epsilon);

    std::cout << "Read Time: " << timer.elapsed() << " seconds" << std::endl;

    PolygonSlicer slicer;

    timer.reset();
    auto layers = slicer.sliceModel(mesh, layerHeight);
    std::cout << "Slice Time: " << timer.elapsed() << " seconds" << std::endl;

    timer.reset();
    rasterizeLayers(mesh, layers);
    std::cout << "Rasterize Time: " << timer.elapsed() << " seconds" << std::endl;

//    minimizeContour(layers[layers.size()-1]->contours[0]);

//    if (true) {
//        print(layers);
//    }
}
