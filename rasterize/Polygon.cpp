//
// Created by erainero on 11/20/18.
//

#include <iostream>
#include <fstream>
#include <cstdint>
#include <cmath>

#include "Polygon.h"
#include "JsonUtils.h"

#if 0


Polygon2i CreatePolygonFromJSON(Json::Value& points) {
    Polygon2i polygon;
    Point2i p;
    Point2i lastp;

    for (int j = 0; j < points.size(); j++) {
        Json::Value point1 = points[j];
        p.x = round(point1[0].asDouble());
        p.y = round(point1[1].asDouble());
        if (j > 0 && !(lastp.x == p.x && lastp.y == p.y)) {
            polygon.pt.push_back(p);
        }
        lastp = p;
    }

    // remove horizontal or vertical redundant points
    for (int j = polygon.pt.size() - 2; j > 0; j--) {
        if (
                (polygon.pt[j - 1].x == polygon.pt[j].x && polygon.pt[j].x == polygon.pt[j + 1].x) ||
                (polygon.pt[j - 1].y == polygon.pt[j].y && polygon.pt[j].y == polygon.pt[j + 1].y)
                )
            polygon.pt.erase(polygon.pt.begin() + j);
    }
    return polygon;
}


#endif
