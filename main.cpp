#include "PolygonSlicer.h"
#include "Timer.h"

int main() {
//    TriangleMesh mesh("../cube_10x10x10.stl", 0.0001);
    TriangleMesh mesh("../bblocky.stl", 0.0001);

    PolygonSlicer slicer;

    Timer timer;
    auto layers = slicer.sliceModel(mesh, 1.0);
    std::cout << "Time: " << timer.elapsed() << std::endl;

    std::cout << "Slices: " << layers.size() << std::endl;

    for (auto layer: layers) {
        std::cout << "Layer: " << layer->z << std::endl;

        for (auto contour: layer->contours) {
            std::cout << "\tContour - external: " << contour.external << " - clockwise: " << contour.clockwise << std::endl;

            for (auto point: contour.points) {
                std::cout << "\t\t" << point << std::endl;
            }
        }
    }
}
