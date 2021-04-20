#ifndef _TREESEPARATION_H_
#define _TREESEPARATION_H_

#include "gurobi_c++.h"
#include "Graph.h"

class TreeSeparation : public GRBCallback
{
public:
    GRBEnv *env;
    GRBVar *y;
    GRBVar *x;
    Graph *G;
    int n;
    int m;
    int r;
    int p;
    int userCuts;
    int lazyCuts;
    int branch_tree_cuts;
    int lp_tree_cuts;
    double userTime;
    double lazyTime;
    TreeSeparation(GRBEnv *env_, GRBVar *yvars, GRBVar *xvars, Graph &graph, int size, int p_);

protected:
    void callback();                            // Path separation callback
    void printSol(double *ybar, double *xbar); // Helper function to print current solution
    void treeCutLP(int root, const double *xbar, const std::vector<std::vector<int>> &tree_adj, const std::vector<int> &num_children);
    void compute_branch_number(const double *xbar, std::vector< std::vector<int> > tree_adj, std::vector<int> tree_degree);
    void compute_tree_weight(const std::vector<std::vector<int>> &tree_adj, std::vector<bool> &discovered, std::vector<int> &s, std::vector<int> &t, int &w, int root);
};

#endif