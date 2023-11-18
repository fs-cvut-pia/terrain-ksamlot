//
// Created by samlo on 18/11/2023.
//

#include "BoatPath.h"

BoatPath::BoatPath(TerrainMap &m, const std::string &nameIn, const Point &startIn, const Point &finishIn) :
        Path(m,nameIn,startIn,finishIn) {}

bool BoatPath::find() {
    return false;
}


