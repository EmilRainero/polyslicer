//
// Created by erainero on 11/20/18.
//

#ifndef RASTERIZE_NODE_H
#define RASTERIZE_NODE_H


class Node {
public:
    int yUpper;
    float xIntersect, dxPerScan;

    Node() : yUpper(-1), xIntersect(-1.0), dxPerScan(0.0) {
        static int idCounter = 0;
        id = idCounter++;
    };

    bool operator< (const Node & b) const
    {
        return (xIntersect < b.xIntersect) ||
                (xIntersect == b.xIntersect && id < b.id);
    }

    bool operator= (const Node & b) const
    {
        return (xIntersect == b.xIntersect && id ==b.id);
    }
private:
    int id;
};


#endif //RASTERIZE_NODE_H
