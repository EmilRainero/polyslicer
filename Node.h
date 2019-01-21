//
// Created by erainero on 11/20/18.
//

#ifndef RASTERIZE_NODE_H
#define RASTERIZE_NODE_H


class Node {
public:
    Node() : yUpper(-1), xIntersect(-1.0), dxPerScan(0.0) {};
    int yUpper;
    float xIntersect, dxPerScan;

    bool operator< (const Node & b) const
    {
        return (xIntersect < b.xIntersect);
    }
};


#endif //RASTERIZE_NODE_H
