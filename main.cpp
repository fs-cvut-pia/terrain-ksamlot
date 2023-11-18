#include "TerrainMap.h"
#include "Path.h"
#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include "PlanePath.h"
#include "BoatPath.h"

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
    const int nx = 256;
    const int ny = 256;

    std::string terrain_filename;

    // Load the terrain map

    if (argc > 1) terrain_filename = argv[1];
    else { std::cout << "No terrain file specified!" << std::endl; return 0; }

    TerrainMap m(nx,ny,terrain_filename);

    // Load the coordinates of the start and end points

    Point start = read_coordinates(argc,argv,2); //pustim terrain.dat 10 10 22 22, takze i_option 2 je prvni 10
    Point finish = read_coordinates(argc,argv,4); //pustim terrain.dat 10 10 22 22, takze i_option 4 je prvni 22

/*    PlanePath planePath(m, start, finish);
    if (planePath.find()) {
        planePath.printStats();
        planePath.saveToFile();
    } else {
        std::cout << "No path found." << std::endl;
    }*/

    std::vector<std::unique_ptr<Path>> paths;
    paths.push_back(std::make_unique<PlanePath>(m,start,finish));

    for (auto& p : paths) {
        std::cout << "Path search: " << p->getName() << std::endl;
        std::cout << "=============" << std::endl;
        p->find();
        p->printStats(); //information of path
        std::cout << "=============" << std::endl;
        p->saveToFile(); //save information of path
    }

    return 0;
}
