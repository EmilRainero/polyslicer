//
// Created by erainero on 11/20/18.
//

#ifndef RASTERIZE_EDGETABLE_H
#define RASTERIZE_EDGETABLE_H

#include <list>
#include <vector>
#include "Point.h"
#include "Polygon.h"
#include "Node.h"
#include "Image.h"
#include "rasterize/Color.h"


class EdgeTable {
public:
    void buildTable(const Polygon2i &);

    int yNext(int, std::vector <Point2i>);

    void makeEdgeRecord(Point2i, Point2i, int);

    static void fillScan(int y, std::list <Node> L, Image& image, const Color& value);

    static void scanFill(Polygon2i P, Image& image, const Color& value);

    void printEdgeTable();

    std::vector <std::list<Node>> Edges;
};



#endif //RASTERIZE_EDGETABLE_H
