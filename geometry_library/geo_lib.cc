#include <algorithm>
#include "geo_lib.h"

void Geo::findViewLimitVertices(Point t, std::vector<Point> pol, int &i1, int &i2) {
    std::vector<BoundaryNode *> clone;
    for (int i = 0; i < pol.size(); i++) {
//        BoundaryNode b = BoundaryNode(pol[i].x_, pol[i].y_);
        BoundaryNode *b = new BoundaryNode();
        b->x_ = pol[i].x_; b->y_ = pol[i].y_;
        b->id_ = i;
        clone.push_back(b);
    }

    BoundaryNode pivot;
    pivot.x_ = t.x_; pivot.y_ = t.y_;
    std::sort(clone.begin(), clone.end(), POLAR_ORDER(pivot));
    i1 = clone[0]->id_;
    i2 = clone[clone.size() - 1]->id_;
}

double Geo::pathLength(std::vector<BasePathPoint> p) {
    double length = 0;
    for (int i = 0; i < p.size() - 1; i++)
    {
        length += G::distance(p[i], p[i+1]);
    }
    return length;
}

//std::vector<CorePolygon> Geo::findObstacles(std::vector<CorePolygon> corePols, Point s, Point t) {
//    // sap xep core polygon theo thu tu xa s dan
//    std::vector<CorePolygon *> clone;
//    for (int i = 0; i < corePols.size(); i++)
//        clone.push_back(&corePols[i]);
//    std::sort(clone.begin(), clone.end(), CORE_DIST_ORDER(s));
//
//    // find obstacles
//    std::vector<CorePolygon> result;
//    for (int i = 0; i < clone.size(); i++) {
//        if (G::segmentPolygonIntersect(s, t, *clone[i]));
//        result.push_back(*clone[i]);
//    }
//
//    return result;
//}
