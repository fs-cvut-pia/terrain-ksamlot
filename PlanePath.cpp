//
// Created by samlo on 18/11/2023.
//

#include "PlanePath.h"
#include <algorithm>
#include <queue>
#include <set>
#include <map>
#include <cmath>

PlanePath::PlanePath(TerrainMap &m, const Point &startIn, const Point &finishIn)
        : Path(m, "Plane", startIn, finishIn) {}

bool PlanePath::find() {
    std::queue<Point> queue;
    std::set<Point> visited;
    std::map<Point, Point> predecessor;
    std::map<Point, double> cost; // Add a map to store the cost to reach each point
    ///
    queue.push(start);
    visited.insert(start);
    cost[start] = 0;

    while (!queue.empty()) {
        Point current = queue.front();
        queue.pop();
        if (current == finish) {
            reconstructPath(predecessor);
            return true;
        }
        for (const auto& neighbor : findNeighbor(current)) {
            double newCost = cost[current] + neighbor.second; // The cost from the current node to the neighbor
            if (!visited.count(neighbor.first) || newCost < cost[neighbor.first]) { // Check if new path to neighbor is shorter
                cost[neighbor.first] = newCost;
                queue.push(neighbor.first);
                visited.insert(neighbor.first);
                predecessor[neighbor.first] = current;
            }
        }
    }
    return false;
}

std::vector<std::pair<Point, double>> PlanePath::findNeighbor(const Point &current) {
    std::vector<std::pair<Point, double>> neighbors;
    for (int j = -1; j < 2; ++j) {
        for (int i = -1; i < 2; ++i) {
            auto neighbor = Point(current.x + i, current.y + j);
            double cost = (j == 0 || i == 0) ? 1.0 : std::sqrt(2.0); // If not diagonal, cost is 1, otherwise sqrt(2)
            if (isValid(neighbor)) {
                neighbors.push_back(std::make_pair(neighbor, cost));
            }
        }
    }
    return neighbors;
}


bool PlanePath::isValid(const Point &referencePoint) {
    return map.validCoords(referencePoint);
}

void PlanePath::reconstructPath(const std::map<Point, Point> &predecessor){
    Point current = finish;
    while (current != start){
        path.push_back(current);
        current = predecessor.at(current);
    }
    path.push_back(current); //pridani startu
    std::reverse(path.begin(), path.end());

}