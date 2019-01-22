//
// Created by erainero on 11/20/18.
//

#include <iostream>

#include "EdgeTable.h"

Node EmptyNode;  // an empty node
std::list <Node> EmptyList;  // an empty list

void insertEdge(std::list <Node> &orderedList, const Node &item) {
    std::list<Node>::iterator curr = orderedList.begin(),
            stop = orderedList.end();
    while ((curr != stop) && ((*curr).xIntersect < item.xIntersect))
        curr++;
    orderedList.insert(curr, item);
}


void writeListInfo(std::list <Node> &L) {
    std::list<Node>::iterator iter;
    for (iter = L.begin(); iter != L.end(); iter++)
        std::cout << "  contents: " << (*iter).yUpper << " "
                  << (*iter).xIntersect << " " << (*iter).dxPerScan << std::endl;
}

void EdgeTable::printEdgeTable() {
    for (int i = 0; i < Edges.size(); i++) {
        if (Edges[i].size() > 1) {
            std::cout << "Scan Line: " << i << " Information" << std::endl;
            writeListInfo(Edges[i]);
        }
    }
}

void buildAEL(std::list <Node> &AEL, std::list <Node> ET) {
    std::list<Node>::iterator iter;

    iter = ET.begin();
    // every Edge table list has a "empty" node at front
    iter++;
    while (iter != ET.end()) {
        insertEdge(AEL, *iter);
        iter++;
    }
}

void EdgeTable::fillScan(int y, std::list <Node> L, Image& image, const Color& value) {
    // want to pull off pairs of x values from adjacent
    // nodes in the list - the y value = scan line
    std::list<Node>::iterator iter1 = L.begin(), iter2;

    int x1, x2;
    while (iter1 != L.end()) {
        iter2 = iter1;
        iter2++;
        x1 = (int) (*iter1).xIntersect;
        x2 = (int) (*iter2).xIntersect;
        image.setPixels(y, x1, x2, value);
        // move on to next pair of nodes
        iter1 = iter2;
        iter1++;
    }
}

void updateAEL(int y, std::list <Node> &L) {   // delete completed edges
    // update the xIntersect field
    std::list<Node>::iterator iter = L.begin();
    while (iter != L.end())
        if (y >= (*iter).yUpper)
            L.erase(iter++);
        else {
            (*iter).xIntersect += (*iter).dxPerScan;
            iter++;
        }
}

void resortAEL(std::list <Node> &L) {
    Node n;
    std::list <Node> L1;
    std::list<Node>::iterator iter = L.begin();
    // create a new list from the old
    // note that the sort command for a list would
    // need us to overload the comparisons operators in the
    // Node class. This is probably just as simple
    while (iter != L.end()) {
        insertEdge(L1, *iter);
        L.erase(iter++);
    }
    L = L1;
}

int EdgeTable::yNext(int k, std::vector <Point2i> p) {
    int j;
    // next subscript in polygon
    if ((k + 1) > (p.size() - 1))
        j = 0;
    else
        j = k + 1;
    //ER
    while (p[k].y == p[j].y)
        if ((j + 1) > (p.size() - 1))
            j = 0;
        else
            j++;
    return (p[j].y);
}


void EdgeTable::makeEdgeRecord(Point2i lower, Point2i upper,
                               int yComp) {
    Node n;

    n.dxPerScan = (float) (upper.x - lower.x) / (upper.y - lower.y);
    n.xIntersect = lower.x;
//    if (upper.y < yComp)
//        n.yUpper = upper.y - 1;
//    else
    n.yUpper = upper.y;
//    std::cout << "EDGE " << n.yUpper << " " << n.xIntersect << " " << n.dxPerScan << std::endl;
    insertEdge(Edges[lower.y], n);
}

void EdgeTable::buildTable(const Polygon2i &Poly) {
    Point2i v1, v2;
    int i, yPrev;

    yPrev = Poly.pt[Poly.pt.size() - 2].y;
    v1.x = Poly.pt[Poly.pt.size() - 1].x;
    v1.y = Poly.pt[Poly.pt.size() - 1].y;
    for (i = 0; i < Poly.pt.size(); i++) {
        v2 = Poly.pt[i];
        if (v1.y != v2.y) { // non horizontal edge
            if (v1.y < v2.y)
                makeEdgeRecord(v1, v2, yNext(i, Poly.pt)); //up edge
            else
                makeEdgeRecord(v2, v1, yPrev); // down edge
        }
        yPrev = v1.y;
        v1 = v2;
    }
}

void EdgeTable::scanFill(Polygon2i P, Image& image, const Color& value) {   // need an edge table and AEL
    EdgeTable EdgeTable;
    std::list <Node> AEL;

//    for (auto pt: P.pt) {
//        std::cout << "  " << pt.x << "," << pt.y << std::endl;
//    }
    EmptyList.push_front(EmptyNode); // and empty list
    // build the edge table - need the window size
    for (int i = 0; i < image.height; i++)
        EdgeTable.Edges.push_back(EmptyList);
    EdgeTable.buildTable(P);
    // if needed - print the table here
//    EdgeTable.printEdgeTable();
    // filling requires building and using AEL
    for (int scanLine = 0; scanLine < image.height; scanLine++) {    // could add output data on table
        buildAEL(AEL, EdgeTable.Edges[scanLine]);
        if (!AEL.empty()) {    // if needed print the table
//            std::cout << "SCANLINE " << scanLine << std::endl;
//            writeListInfo(AEL);
//            std::cout << "Scan line: " << scanLine << " " << AEL.size() << std::endl;



//            writeListInfo(AEL);
            updateAEL(scanLine, AEL);

//            std::cout << "updated" << std::endl;

            fillScan(scanLine, AEL, image, value);
//            std::cout << "  filled: " << AEL.size() << std::endl;

            resortAEL(AEL);
//            std::cout << "  resorted: " << AEL.size() << std::endl;

//            std::cout << "Done" << std::endl;
        }
    }
}
