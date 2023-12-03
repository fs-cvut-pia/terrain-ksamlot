//
// Created by samlo on 21/11/2023.
//

#include "CarPath.h"
#include <algorithm>
#include <queue>
#include <set>


CarPath::CarPath(TerrainMap &m, const Point &startIn, const Point &finishIn)
        : Path(m,"Car",startIn,finishIn) {}

bool CarPath::find() {
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

std::vector<Point> CarPath::findNeighbor(const Point &current) {
    std::vector<Point> neighbors;
    for (int j = -1; j < 2; ++j) {
        for (int i = -1; i < 2; ++i) {
            auto neighbor = Point(current.x + i, current.y + j);
            if (isValid(neighbor, current)) {
                //std::cout<<"For point: ["<<current.x << "," << current.y<<"] neighbor: [" << neighbor.x << "," << neighbor.y<<"]" << "altitude: "<< map.alt(neighbor)<<std::endl;
                neighbors.push_back(neighbor);
            }
        }
    }
    return neighbors;
}

bool CarPath::isValid(const Point &referencePoint, const Point &currentPoint) {
    if (map.validCoords(referencePoint)){
        if(map.alt(referencePoint)>0 || referencePoint == finish) {
            // Ensure that the altitude difference between reference and current point is <= 0.6%
            if(abs((map.alt(referencePoint) - map.alt(currentPoint)) / map.alt(currentPoint)) <= 0.006 || referencePoint == finish){
                return true;
            }
        }
    }
    return false;
}


/*std::vector<Point> CarPath::findNeighbor(const Point &current) {
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

bool CarPath::isValid(const Point &referencePoint) {
    if (map.validCoords(referencePoint)){
        if(map.alt(referencePoint)>0 || referencePoint == finish)
            return true;
    }
    return false;
}*/

void CarPath::reconstructPath(const std::map<Point, Point>& predecessor){
    Point current = finish;
    while (current != start){
        path.push_back(current);
        current = predecessor.at(current);
    }
    path.push_back(current); //pridani startu
    std::reverse(path.begin(), path.end());
}