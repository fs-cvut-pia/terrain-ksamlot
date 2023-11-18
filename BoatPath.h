//
// Created by samlo on 18/11/2023.
//

#ifndef TERRAIN_KSAMLOT_BOATPATH_H
#define TERRAIN_KSAMLOT_BOATPATH_H

#include "Path.h"

class BoatPath : public Path{
public:
    BoatPath(TerrainMap &m, const std::string &nameIn, const Point &startIn, const Point &finishIn);

    bool find() override;
};

#endif //TERRAIN_KSAMLOT_BOATPATH_H
