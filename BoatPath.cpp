//
// Created by samlo on 18/11/2023.
//

#include "BoatPath.h"
#include <iostream>
#include <algorithm>
#include <queue>
#include <set>
#include <map>

BoatPath::BoatPath(TerrainMap &m, const Point &startIn, const Point &finishIn)
        : Path(m,"Boat",startIn,finishIn) {}

bool BoatPath::find() {
    std::queue<Point> queue;
    std::set<Point> visited;
    std::map<Point, Point> predecessor;
    ///
    queue.push(start);
    visited.insert(start);

    while (!queue.empty()) {
        Point current = queue.front();
        queue.pop();
        for (const auto neighbor: findNeighbor(current)) {
            if (visited.find(neighbor) == visited.end()){
                queue.push(neighbor);
                visited.insert(neighbor);
                predecessor[neighbor] = current; //mapa(predicessor) na pozici klice[] = hodnota
                if(neighbor == finish){
                    reconstructPath(predecessor);
                    return true;
                }
            }
        }
    }
    return false;
}

std::vector<Point> BoatPath::findNeighbor(const Point &current) {
    std::vector<Point> neighbors;
    for (int j = -1; j < 2; ++j) {
        for (int i = -1; i < 2; ++i) {
            auto neighbor = Point(current.x + i, current.y + j);
            if (isValid(neighbor)) {
                //std::cout<<"Pro bod: ["<<current.x << "," << current.y<<"] soused: [" << neighbor.x << "," << neighbor.y<<"]" << "výška: "<< map.alt(neighbor)<<std::endl;
                neighbors.push_back(neighbor);
            }
        }
    }
    return neighbors;
}

bool BoatPath::isValid(const Point &referencePoint) {
    if (map.validCoords(referencePoint)){
        if(map.alt(referencePoint)<0 || referencePoint == finish)
        return true;
    }
    return false;
}
//bool BoatPath::isValid(const Point &referencePoint) {
//    if (referencePoint.x >= map.nx || referencePoint.x < 0 ||
//        referencePoint.y >= map.ny || referencePoint.y < 0) {
//        return false;
//    }
//    try {
//        if(map.alt(referencePoint)>0)
//            return false;
//    }
//    catch (const std::exception& e) {
//        std::cout << "Exception caught: " << e.what() << std::endl;
//    }
//    return true;
//}

void BoatPath::reconstructPath(const std::map<Point, Point>& predecessor){
    Point current = finish;
    while (current != start){
        path.push_back(current);
        current = predecessor.at(current);
    }
    path.push_back(current); //pridani startu
    std::reverse(path.begin(), path.end());
}
