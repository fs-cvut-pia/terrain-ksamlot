//
// Created by samlo on 18/11/2023.
//

#include "PlanePath.h"
#include <algorithm>
#include <queue>
#include <set>
#include <map>
#include <algorithm>


PlanePath::PlanePath(TerrainMap &m, const Point &startIn, const Point &finishIn)
        : Path(m, "Letadlo", startIn, finishIn) {}

bool PlanePath::find() {
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

std::vector<Point> PlanePath::findNeighbor(const Point &current) {
    std::vector<Point> neighbors;
    for (int j = -1; j < 2; ++j) {
        for (int i = -1; i < 2; ++i) {
            auto neighbor = Point(current.x + i, current.y + j);
            if (isValid(neighbor)) {
                neighbors.push_back(neighbor);
            }
        }
    }
    return neighbors;
}

bool PlanePath::isValid(const Point &referencePoint) {
    if (referencePoint.x > map.nx || referencePoint.x < 0 || referencePoint.y > map.ny || referencePoint.y < 0) {
        return false;
    }
    return true;
}

void PlanePath::reconstructPath(const std::map<Point, Point>& predecessor){
    Point current = finish;
    while (current != start){
        path.push_back(current);
        current = predecessor.at(current);
    }
    path.push_back(current); //pridani startu
    std::reverse(path.begin(), path.end());

}