#ifndef _PATHSEPARATION_H_
#define _PATHSEPARATION_H_

#include "gurobi_c++.h"
#include "Graph.h"

class PathSeparation : public GRBCallback
{
public:
    GRBVar **y;
    GRBVar *x;
    Graph *G;
    int n;
    int m;
    int r;
    int userCuts;
    int lazyCuts;
    double userTime;
    double lazyTime;
    PathSeparation(GRBVar **yvars, GRBVar *xvars, Graph &graph, int size);

protected:
    void callback();                            // Path separation callback
    void printSol(double **ybar, double *xbar); // Helper function to print current solution
};

#endif