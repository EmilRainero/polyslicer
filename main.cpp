#include "PolygonSlicer.h"
#include "Timer.h"

int main() {
    TriangleMesh mesh("../cube_10x10x10.stl", 0.0001);

    PolygonSlicer slicer;

    Timer timer;
    slicer.sliceModel(mesh, 0.1);
    std::cout << "Time: " << timer.elapsed() << std::endl;
}
