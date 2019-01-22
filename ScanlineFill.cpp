#include "ScanlineFill.h"
#include <iostream>
#include <mutex>

#include "EdgeTable.h"

std::mutex scanfilllock;


void ScanlineFill::insertEdge(std::set <Node> &orderedList, const Node &item) {
//    std::list<Node>::iterator curr = orderedList.begin(),
//            stop = orderedList.end();
//    while ((curr != stop) && ((*curr).xIntersect < item.xIntersect))
//        curr++;
//    orderedList.insert(curr, item);

//    for (auto n: orderedList) {
//        std::cout << n.xIntersect << " ";
//    }
//    std::cout << "insert " << item.xIntersect << std::endl;
    orderedList.insert(item);
}


void ScanlineFill::writeListInfo(std::set <Node> &L) {
    for (auto node: L) {
        std::cout
                << "  contents: "
                   << node.yUpper << " "
                     << node.xIntersect << " "
                << node.dxPerScan
                  << std::endl;
    }
}

void ScanlineFill::printEdgeTable() {
    for (int i = 0; i < Edges.size(); i++) {
        if (Edges[i].size() > 1) {
            std::cout << "Scan Line: " << i << " Information" << std::endl;
            writeListInfo(Edges[i]);
        }
    }
}

void ScanlineFill::buildAEL(std::set <Node> &AEL, std::set <Node> ET) {
//    std::list<Node>::iterator iter;
//
//    iter = ET.begin();
//    // every Edge table list has a "empty" node at front
//    iter++;
//    while (iter != ET.end()) {
//        insertEdge(AEL, *iter);
//        iter++;
//    }

//    std::cout << "buildAEL before " << AEL.size() << "  ET " << ET.size()  << std::endl;
    for (auto node: ET) {
        insertEdge(AEL, node);
    }
//    std::cout << "buildAEL after " << AEL.size() << std::endl;
}

void ScanlineFill::fillScan(int y, std::set <Node> L, Image& image, const Color& value) {
    // want to pull off pairs of x values from adjacent
    // nodes in the list - the y value = scan line
    std::set<Node>::iterator iter1 = L.begin(), iter2;

    if (L.size() % 2 != 0) {
        std::cout << "fillScan # nodes " << L.size() << " layer: " << layerNumber << " y: " << y << std::endl;
        return;
    }

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

void ScanlineFill::updateAEL(int y, std::set <Node> &L) {   // delete completed edges
    // update the xIntersect field
//    std::set<Node>::iterator iter = L.begin();
//    while (iter != L.end())
//        if (y >= (*iter).yUpper)
//            L.erase(iter++);
//        else {
//            (*iter).xIntersect += (*iter).dxPerScan;
//            iter++;
//        }


    std::list<std::set<Node>::iterator> nodesToDelete;

//    std::cout << "updateAEL before " << L.size() << " items" << std::endl;

    for (std::set<Node>::iterator iter = L.begin(); iter != L.end(); ++iter) {
        if (y == (*iter).yUpper)
            nodesToDelete.push_back(iter);
        else {
            Node &writableNode = const_cast<Node&>(*iter);
            writableNode.xIntersect += writableNode.dxPerScan;
        }
    }
    for (auto nodeIter: nodesToDelete) {
//        std::cout << "erase " << node.yUpper << std::endl;
        L.erase(nodeIter);
    }

//    std::cout << "updateAEL after " << L.size() << " items" << std::endl;

//    for (Node& node: L) {
//        if (y >= node.yUpper) {
//
//        } else {
//            node.xIntersect = node.xIntersect + node.dxPerScan;
//        }
//    }
//    for (auto& it = L.begin(); it != L.end(); ) {
//        if (y >= (*it).yUpper)
//            L.erase(it++);
//        else {
//            (*it).xIntersect += (*it).dxPerScan;
//            it++;
//        }
//    }
}

void ScanlineFill::resortAEL(std::set <Node> &L) {
//    Node n;
//    std::list <Node> L1;
//    std::list<Node>::iterator iter = L.begin();
//    // create a new list from the old
//    // note that the sort command for a list would
//    // need us to overload the comparisons operators in the
//    // Node class. This is probably just as simple
//    while (iter != L.end()) {
//        insertEdge(L1, *iter);
//        L.erase(iter++);
//    }
//    L = L1;
}

int ScanlineFill::yNext(int k, std::vector <Point2i> p) {
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


void ScanlineFill::makeEdgeRecord(Point2i lower, Point2i upper,
                               int yComp) {
    Node n;

    n.dxPerScan = (float) (upper.x - lower.x) / (upper.y - lower.y);
    n.xIntersect = lower.x;
//    if (upper.y < yComp)
//        n.yUpper = upper.y - 1;
//    else
    n.yUpper = upper.y;
//    std::cout << "EDGE " << n.yUpper << " " << n.xIntersect << " " << n.dxPerScan << std::endl;

//    insertEdge(Edges[lower.y], n);

    Edges[lower.y].insert(n);

}

void ScanlineFill::buildTable(const Polygon2i &Poly) {
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


void ScanlineFill::scanFill(Polygon2i P, Image& image, const Color& value, int layerNumber) {
    ScanlineFill edgeTable;
    std::set <Node> AEL;

    this->layerNumber = layerNumber;
    std::set <Node> EmptyList;  // an empty list

    // print polygon
//    for (auto pt: P.pt) {
//        std::cout << "  " << pt.x << "," << pt.y << std::endl;
//    }
    for (int i = 0; i < image.height; i++)
        edgeTable.Edges.push_back(EmptyList);

    scanfilllock.lock();
    edgeTable.buildTable(P);
    scanfilllock.unlock();

//    edgeTable.printEdgeTable();
    // filling requires building and using AEL
    for (int scanLine = 0; scanLine < image.height; scanLine++) {    // could add output data on table
//        std::cout << "Scanline: " << scanLine << " edgetable size: " << edgeTable.Edges[scanLine].size() << std::endl;
        buildAEL(AEL, edgeTable.Edges[scanLine]);
        if (!AEL.empty()) {    // if needed print the table
            updateAEL(scanLine, AEL);
            fillScan(scanLine, AEL, image, value);
//            resortAEL(AEL);
        }
    }
}
