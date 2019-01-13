#include <cmath>
#include "PolygonSlicer.h"
#include "Timer.h"


void print(std::vector<Layer *> layers) {
    for (auto layer: layers) {
        std::cout << "Layer: " << layer->z << std::endl;

        for (auto contour: layer->contours) {
            std::cout << "\tContour - external: " << contour.external << " - clockwise: " << contour.clockwise
                      << std::endl;

            for (auto point: contour.points) {
                std::cout << "\t\t" << point << std::endl;
            }
        }
    }
}

void rescaleLayer(Layer* layer, float tx, float ty, float scale) {
    for (Contour& contour: layer->contours) {
        for (Vector3& point: contour.points) {
            point.x = (point.x - tx) * scale;
            point.y = (point.y - ty) * scale;
        }
    }
}

void scaleLayers(TriangleMesh& mesh, std::vector<Layer *> layers) {
    auto minVertex = mesh.getBottomLeftVertex();
    auto maxVertex = mesh.getTopRightVertex();

    float xSize = maxVertex.x - minVertex.x;
    float ySize = maxVertex.y - minVertex.y;
    float maxScale = std::max(xSize, ySize);

    std::cout << xSize << " " << ySize << " " << maxScale << std::endl;

    float arraySize{1024};
    float scale = (arraySize - 1) / maxScale;
    for (auto layer: layers) {
        rescaleLayer(layer, minVertex.x, minVertex.y, scale);
    }
}

int main() {
    double epsilon{0.0001};
    double layerHeight{0.1};

    TriangleMesh mesh("../cube_10x10x10.stl", epsilon);
//    TriangleMesh mesh("../bblocky.stl", epsilon);
//    TriangleMesh mesh("../pisa.stl", epsilon);

    PolygonSlicer slicer;

    Timer timer;
    auto layers = slicer.sliceModel(mesh, layerHeight);
    std::cout << "Time: " << timer.elapsed() << " seconds" << std::endl;
    std::cout << "Slices: " << layers.size() << std::endl;

    scaleLayers(mesh, layers);

    if (true) {
        print(layers);
    }
}
