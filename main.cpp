#include <iostream>
#include <vector>
#include "geometry_library/geo_lib.h"
#include "graph/o_graph.h"

using namespace std;

int main() {
    std::vector<CorePolygon> core_set;
    std::vector<Point> vertices;
    CorePolygon tmp = CorePolygon();

    vertices.push_back(Point(4.09, 57.14));
    tmp.vertices = vertices;
    tmp.id = 0;

    core_set.push_back(tmp);

    vertices.clear();
    vertices.push_back(Point(13.22, 44.25));
    vertices.push_back(Point(4.81, 29.22));
    vertices.push_back(Point(17.69, 7.75));
    vertices.push_back(Point(56.71, 26.18));
    vertices.push_back(Point(36.66, 50.34));
    tmp.vertices = vertices;
    tmp.id = 1;

    core_set.push_back(tmp);

    vertices.clear();
    vertices.push_back(Point(51.34, -4.07));
    vertices.push_back(Point(69.24, 20.81));
    vertices.push_back(Point(98.94, 11.32));
    vertices.push_back(Point(88.21, -25.36));
    tmp.vertices = vertices;
    tmp.id = 2;

    core_set.push_back(tmp);

    vertices.clear();
    vertices.push_back(Point(121.31, -22.5));
    tmp.vertices = vertices;
    tmp.id = 3;

    core_set.push_back(tmp);

    OGraph *g = new OGraph(core_set);
    std::vector<std::vector<BasePathPoint>> all_paths_ = g->basePaths(0.5, 8);

    for (int i = 0; i < all_paths_.size(); i++)
    {
        for (int j = 0; j < all_paths_[i].size(); j++)
            printf("(%d, %g, %g)\t", all_paths_[i][j].hole_id_, all_paths_[i][j].x_, all_paths_[i][j].y_);
        printf("\n");
    }

    printf("%g\n", g->shortestPathLength());

    return 0;

}