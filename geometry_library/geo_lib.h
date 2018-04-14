#pragma once
#include "../geomathhelper/geo_math_helper.h"

class Geo;

struct BoundaryNode : Point {
    int id_;
    bool is_core_polygon_boundary;

    // comparison is done first on y coordinate and then on x coordinate
    bool operator<(BoundaryNode n2) {
        return y_ == n2.y_ ? x_ < n2.x_ : y_ < n2.y_;
    }

    BoundaryNode(double x, double y) {
        Point(x, y);
    }
    BoundaryNode() {}
};

struct BasePathPoint: Point {
    int hole_id_;
};

// Hole struct
struct HoleSt {
    std::vector<BoundaryNode> node_vector_; // node vector of hole
    int id_;   // id of hole
};

// core polygon
struct CorePolygon {
    std::vector<Point> vertices;
    int id;
};

struct Tangent {
    Point start;
    Point end;
    double length;
    int start_hole_id_;
    int end_hole_id_;

    Tangent(Point a, Point b) {
        start = a;
        end = b;
        length = G::distance(a, b);
    }

    Tangent(Point a, Point b, int start_id_, int end_id_) {
        start = a;
        end = b;
        length = G::distance(a, b);
        start_hole_id_ = start_id_;
        end_hole_id_ = end_id_;
    }
};

struct COORDINATE_ORDER {
    bool operator()(const BoundaryNode *a, const BoundaryNode *b) const {
        return a->y_ == b->y_ ? a->x_ < b->x_ : a->y_ < b->y_;
    }
    bool operator()(const node *a, const node *b) const {
        return a->y_ == b->y_ ? a->x_ < b->x_ : a->y_ < b->y_;
    }

};

// used for sorting points according to polar order w.r.t the pivot
struct POLAR_ORDER {
    POLAR_ORDER(struct BoundaryNode p) { this->pivot = p; }

    bool operator()(const BoundaryNode *a, const BoundaryNode *b) const {
        int order = G::orientation(pivot, *a, *b);
        if (order == 0)
            return G::distance(pivot, *a) < G::distance(pivot, *b);
        return (order == 2);
    }

    struct BoundaryNode pivot;
};

// used for sorting corepolygons according to distance from source
struct CORE_DIST_ORDER {
    CORE_DIST_ORDER(Point p) { this->pivot = p; }
    bool operator()(const CorePolygon *a, const CorePolygon *b) const {
        return G::distance(a->vertices[0], pivot) < G::distance(b->vertices[0], pivot) ? 1 : 0;
    }
    Point pivot;
};

class Geo {
public:
    static void findViewLimitVertices(Point, std::vector<Point>, int &i1, int &i2);
    static double pathLength(std::vector<BasePathPoint> );
//    static std::vector<CorePolygon> findObstacles(std::vector<CorePolygon>, Point, Point);
};