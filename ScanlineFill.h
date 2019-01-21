//
// Created by erainero on 1/21/19.
//

#ifndef POLYSLICER_SCANLINEFILL_H
#define POLYSLICER_SCANLINEFILL_H

#include <list>
#include <set>
#include <vector>
#include "Point.h"
#include "Polygon.h"
#include "Node.h"
#include "Image.h"
#include "Color.h"


class ScanlineFill {
public:
    void buildTable(const Polygon2i &);

    int yNext(int, std::vector <Point2i>);

    void makeEdgeRecord(Point2i, Point2i, int);

    static void fillScan(int y, std::set <Node> L, Image& image, const Color& value);

    static void scanFill(Polygon2i P, Image& image, const Color& value);

    static void insertEdge(std::set <Node> &orderedList, const Node &item);

    static void writeListInfo(std::set <Node> &L);

    static void buildAEL(std::set <Node> &AEL, std::set <Node> ET);

    static void updateAEL(int y, std::set <Node> &L);

    static void resortAEL(std::set <Node> &L);

    void printEdgeTable();

    std::vector <std::set<Node>> Edges;
};


#endif //POLYSLICER_SCANLINEFILL_H