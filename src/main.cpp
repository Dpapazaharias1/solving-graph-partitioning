#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include "Graph.hpp"
#include "Formulation.hpp"

int main(int argc, char *argv[])
{
    if(argc < 5) {
        std::cout << "Not enough arguments. Usage as follows" << std::endl;
        std::cout << "argv[1]: filetype - input file type edge list representation (-e) or adjacency list representation (-a)" << std::endl;
        std::cout << "argv[2]: weighted - unweighted (-u), edge weighted (-ew), randomly generated weights (-w)" << std::endl;
        std::cout << "argv[3]: filename - the instance we wish to solve" << std::endl;
        std::cout << "argv[4]: formulation - Options : -TRI, -TRILP, -FLOW, -FLOW+, -PATH, -PATH+, -TCFLP, -TDPLP, -TDPBLP" << std::endl;
        std::cout << "Depending on formulation selected, addition arguments may be required." << std::endl;
    }
    else {
        const char *filetype = argv[1];
        const char *weighted = argv[2];
        const char *filename = argv[3];
        const char *formulation = argv[4];

        Graph G = Graph(filetype, filename, weighted);

        if(strcmp(formulation, "-TRI") == 0) {
            std::cout << "Starting triangle integer programming formulation." << std::endl;
            triangle(G, true);
        }
        if(strcmp(formulation, "-TRILP") == 0) {
            std::cout << "Starting triangle linear programming formulation." << std::endl;
            triangle(G, false);
        }
        if(strcmp(formulation, "-FLOW") == 0) {
            std::cout << "Starting flow integer programming formulation." << std::endl;
            flow(G, 0, 0, "-none");
        }
        if(strcmp(formulation, "-FLOW+") == 0) {
            std::cout << "Starting flow with cuts integer programming formulation." << std::endl;
            if (argc < 8) {
                std::cout << "Missing arguments for FLOW+ formulation" << std::endl;
                std::cout << "argv[5]: p - probability of generating a user cut at a node" << std::endl;
                std::cout << "argv[6]: r_pct - For knap and tdp cuts; Target weight of the trees to generate cuts" << std::endl;
                std::cout << "argv[7]: cut_type - Fractional separation routine. Options: -prim, -flow, -tdp" << std::endl;
            }
            else {
                int p = std::stoi(argv[5]);
                float r_pct = std::stof(argv[6]);
                const char* cut_type = argv[7];
                flow(G, p, r_pct, cut_type);
            }
        }
        if(strcmp(formulation, "-PATH") == 0) {
            std::cout << "Starting path integer programming formulation." << std::endl;
            path(G, 0, 0, "-path");
        }

        if(strcmp(formulation, "-PATH+") == 0) {
            std::cout << "Starting path with cuts integer programming formulation." << std::endl;
            if (argc < 8) {
                std::cout << "Missing arguments for FLOW+ formulation" << std::endl;
                std::cout << "argv[5]: p - probability of generating a user cut at a node" << std::endl;
                std::cout << "argv[6]: r_pct - For knap and tdp cuts; Target weight of the trees to generate cuts" << std::endl;
                std::cout << "argv[7]: cut_type - Fractional separation routine. Options: -prim, -flow, -tdp" << std::endl;
            }
            else {
                int p = std::stoi(argv[5]);
                float r_pct = std::stof(argv[6]);
                const char* cut_type = argv[7];
                path(G, p, r_pct, cut_type);
            }
        }

        if(strcmp(formulation, "-TCFLP") == 0) {
            std::cout << "Starting tree cover linear programming formulation." << std::endl;
            if (G.m != G.n - 1) { std::cout << "Graph is not a tree." << std::endl;}
            else{ tree_cover_formulation(G);}
        }

        if(strcmp(formulation, "-TCFIP") == 0) {
            std::cout << "Starting tree cover integer programming formulation." << std::endl;
            tree_cover_ip(G, "-tcf");
        }

        if(strcmp(formulation, "-BLOCKIP") == 0) {
            std::cout << "Starting block decomposition integer programming formulation." << std::endl;
            tree_cover_ip(G, "-block");
        }

        if(strcmp(formulation, "-TDPLP") == 0) {
            std::cout << "Starting tree dynamic programming LP formulation." << std::endl;
            if (G.m != G.n - 1) { std::cout << "Graph is not a tree." << std::endl;}
            else { tree_lp(G); }
        }
        if(strcmp(formulation, "-TDPBLP") == 0) {
            std::cout << "Starting tree dynamic programming with Benders cuts LP formulation." << std::endl;
            if (G.m != G.n - 1) { std::cout << "Graph is not a tree." << std::endl;}
            else { dynamic_program_lp(G); }
        }
        if(strcmp(formulation, "-FLOWLP") == 0){
            std::cout << "Starting tree dynamic programming with Benders cuts LP formulation." << std::endl;
            if (argc < 7) {
                std::cout << "Missing arguments for FLOW+ formulation" << std::endl;
                std::cout << "argv[5]: p - probability of generating a user cut at a node" << std::endl;
                std::cout << "argv[6]: r_pct - For knap and tdp cuts; Target weight of the trees to generate cuts" << std::endl;
            }
            else {
                int p = std::stoi(argv[5]);
                float r_pct = std::stof(argv[6]); 
                flow_relax(G, p, r_pct); 
            }
        }
    }
}