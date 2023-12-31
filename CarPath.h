//
// Created by samlo on 21/11/2023.
//

#ifndef PATH_PLOT_PY_CARPATH_H
#define PATH_PLOT_PY_CARPATH_H

#include "Path.h"
#include <map>

class CarPath : public Path{
public:
    CarPath(TerrainMap &m, const Point &startIn, const Point &finishIn);

    bool find() override;

private:
    bool isValid(const Point &referencePoint, const Point &current);

    std::vector<std::pair<Point, double>> findNeighbor(const Point &current); // Change return type here

    void reconstructPath(const std::map<Point, Point> &predecessor);

/*
    bool isValid(const Point& referencePoint, const Point &currentPoint);

    std::vector<Point> findNeighbor(const Point &current);

    void reconstructPath(const std::map<Point, Point> &predecessor);
    */

};

#endif //PATH_PLOT_PY_CARPATH_H
