//
// Created by samlo on 18/11/2023.
//

#ifndef TERRAIN_KSAMLOT_BOATPATH_H
#define TERRAIN_KSAMLOT_BOATPATH_H

#include "Path.h"
#include <map>

class BoatPath : public Path{
public:
    BoatPath(TerrainMap &m, const Point &startIn, const Point &finishIn);

    bool find() override;

private:
    bool isValid(const Point& referencePoint);

    std::vector<std::pair<Point, double>> findNeighbor(const Point &current); // Change return type here

    void reconstructPath(const std::map<Point, Point> &predecessor);

};

#endif //TERRAIN_KSAMLOT_BOATPATH_H
