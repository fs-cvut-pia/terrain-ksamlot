#include "TerrainMap.h"
#include "Path.h"
#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include "PlanePath.h"
#include "BoatPath.h"
#include "CarPath.h"

// Include files of your path classes will need to be added here

Point read_coordinates(int argc, char *argv[], int i_option) {
    //argc - argument count, kolik jich mam +1
    //argv - jednotlivy argumenty, 0 je cesta(nazev)
    Point p;

    if (argc > i_option+1) {
        p.x = std::atoi(argv[i_option]); //pustim terrain.dat 10 10 22 22, takze i_option 2 je prvni 10
        p.y = std::atoi(argv[i_option + 1]); //tohle je i_option +1 takze 3 takze 10
    }
    else throw std::runtime_error("Coordinates incorrectly specified!");

    return p;
}

int main(int argc, char *argv[]) {
    try {
        if (argc <= 4) {
            std::cerr << "Insufficient arguments!" << std::endl;
            return 1;
        }

        const int nx = 256;
        const int ny = 256;

        std::string terrain_filename = argv[1];
        TerrainMap m(nx, ny, terrain_filename);

        Point start = read_coordinates(argc, argv, 2);
        Point finish = read_coordinates(argc, argv, 4);

        // uniqu
        std::vector<std::shared_ptr<Path>> paths = {
                std::make_shared<PlanePath>(m, start, finish),
                std::make_shared<BoatPath>(m, start, finish),
                std::make_shared<CarPath>(m, start, finish)

        };
//        paths.push_back(std::make_unique<BoatPath>(m, start, finish));
//        paths.push_back(std::make_unique<PlanePath>(m, start, finish));

        for (auto& p : paths) {
            std::cout << "Path search: " << p->getName() << std::endl;
            std::cout << "=============" << std::endl;
            p->find();
            p->printStats();
            std::cout << "=============" << std::endl;
            p->saveToFile();
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
