#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

#include "Graph.h"
#include "Formulation.h"


void printResults(char *filename, std::vector<double> solution)
{
    std::cout << filename << " ";
    for(int s = 0; s < 11; s++)
    {
        std::cout << solution[s] << " ";
    }
        std::cout << solution[11] << "\n";
}

void logResults(std::string filepath, std::vector<double> solution)
{
    std::ofstream results;
    std::ifstream ifile(filepath);
    if (ifile)
    {
        // std::cout << "File exists \n";
        results.open(filepath, std::ofstream::app);
        for (int s = 0; s < 11; s++)
        {
            results << solution[s] << ",";
        }
        results << solution[11] << "\n";
        results.close();
    }
    else
    {
        // std::cout << "File does not exist \n";
        results.open(filepath, std::ofstream::app);
        results << "n, m, r, Lazy Cuts, User Cuts, Runtime, Variables, Binaries, Constraints, Nodes, Objective, Gap\n";
        for (int s = 0; s < 11; s++)
        {
            results << solution[s] << ",";
        }
        results << solution[11] << "\n";
        results.close();
    }
}

int main(int argc, char *argv[])
{
    std::cout << "Number of arguments: " << argc << std::endl;
    if (argc == 5)
    {
        const char *filetype = argv[1];
        const char *filename = argv[2];
        const char *weighted = argv[3];
        const char *formulation = argv[4];
        std::cout << "Constructing Graph..." << std::endl;
        Graph G(filetype, filename, weighted);
        if(strcmp(formulation, "-tri") == 0)
        {
            std::cout << "Starting triangle formulation..." << std::endl;
            triangle(G);
        }

        if(strcmp(formulation, "-flow") == 0)
        {
            std::cout << "Starting flow formulation..." << std::endl;
            flow(G);
        }
        if(strcmp(formulation, "-path") == 0)
        {
            std::cout << "Starting path formulation..." << std::endl;
            path(G);
        }

    }
    else
    {
        std::cout << "Error: Missing arguments" << std::endl;
    }
    return 0;
}