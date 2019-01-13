#include "PolygonSlicer.h"
#include "Timer.h"

int main() {
    double epsilon{0.0001};
    double layerHeight{0.1};

//    TriangleMesh mesh("../cube_10x10x10.stl", epsilon);
//    TriangleMesh mesh("../bblocky.stl", epsilon);
    TriangleMesh mesh("../pisa.stl", epsilon);

    PolygonSlicer slicer;

    Timer timer;
    auto layers = slicer.sliceModel(mesh, layerHeight);
    std::cout << "Time: " << timer.elapsed() << " seconds" << std::endl;
    std::cout << "Slices: " << layers.size() << std::endl;

    if (false) {
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
}
