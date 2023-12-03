//
// Created by samlo on 18/11/2023.
//

#ifndef TERRAIN_KSAMLOT_PLANEPATH_H
#define TERRAIN_KSAMLOT_PLANEPATH_H

#include <map>
#include "Path.h"


class PlanePath : public Path{
public:
    PlanePath(TerrainMap &m, const Point &startIn, const Point &finishIn);

    bool find() override;
private:

    bool isValid(const Point& referencePoint);

    std::vector<std::pair<Point, double>> findNeighbor(const Point &current); // Change return type here

    void reconstructPath(const std::map<Point, Point> &predecessor);

};

#endif //TERRAIN_KSAMLOT_PLANEPATH_H

/* Prohledavani do sirky nezohlednujici rozdil ve vzdalenosti pri pohybu sikmo
#ifndef TERRAIN_KSAMLOT_PLANEPATH_H
#define TERRAIN_KSAMLOT_PLANEPATH_H

#include <map>
#include "Path.h"


class PlanePath : public Path{
public:
    PlanePath(TerrainMap &m, const Point &startIn, const Point &finishIn);

    bool find() override;
private:

    bool isValid(const Point& referencePoint);

    std::vector<Point> findNeighbor(const Point &current);

    void reconstructPath(const std::map<Point, Point> &predecessor);

};

#endif //TERRAIN_KSAMLOT_PLANEPATH_H
*/
