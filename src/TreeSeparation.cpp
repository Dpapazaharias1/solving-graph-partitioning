#include <algorithm>
#include <limits.h>
#include <random>
#include <chrono>
#include <queue>
#include "TreeSeparation.h"

TreeSeparation::TreeSeparation(GRBEnv *env_, GRBVar *yvars, GRBVar *xvars, Graph &graph, int size, int p_)
{
    env = env_;
    y = yvars;
    x = xvars;
    r = size;
    n = graph.n;
    m = graph.m;
    p = p_;
    G = &graph;
    userCuts = 0;
    lazyCuts = 0;
    branch_tree_cuts = 0;
    lp_tree_cuts = 0;
    userTime = 0.0;
    lazyTime = 0.0;
}

void TreeSeparation::callback()
{
    try {
        if (where == GRB_CB_MIPNODE) {
            if (getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL){

                std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
                // -------- Retrieve Solution --------

                int rand_int = rand() % 100;

                double *xbar = new double[n * n];
                double *ybar = new double[2 * m];
                xbar = getNodeRel(x, n * n);
                ybar = getNodeRel(y, 2 * m);

                GRBLinExpr expr;
                int *parent = new int[n];
                double *cost = new double[n];
                int j;

                std::vector<double> CostTo(2 * m, 0.0);
                for (int e = 0; e < 2 * m; e++)
                {   
                    CostTo[e] = ybar[e];
                    if (ybar[e] < 0) { CostTo[e] = 0.0; } // Some approximations of 0 are negative
                }
                /*
                for(int i =0 ; i < n; i++)
                {
                    for(int e = G->EdgesBegin[i]; e < G->EdgesBegin[i] + G->Degree[i]; e++)
                    {
                        std::cout << "Edge (" << i << "," << G->EdgeTo[e] << "): " << ybar[e] << std::endl;
                    }
                }
                */
                if(p > rand_int)
                {
                    std::vector<int> treeNodes;
                    srand((unsigned)time(0));
                    int s = rand() % n;
                    double treeWeight = 0.0;
                    int node_weight = 0;
                    treeNodes = G->Prim(CostTo, s, r, parent, cost);
                    
                    expr = 0;
                    for (int i : treeNodes)
                    {
                        node_weight += G->Weight[i];
                        if (parent[i] != -1)
                        {
                            //std::cout << "(" << i << ",";
                            j = parent[i];
                            treeWeight += ybar[j];
                            expr += y[j];
                            //std::cout << G->EdgeTo[j] << "), ";
                        }
                    }
                    
                    //std::cout << "Node weight : " << node_weight << ">" << r << std::endl;
                    //std::cout << "sum x : " << treeWeight << std::endl;

                    
                    if (treeWeight < 1.0 && node_weight > r)
                    {
                        addCut(expr, GRB_GREATER_EQUAL, 1.0);
                        userCuts++;
                    }
                }
                // -------- Memory Dellocation --------
               
                delete[] ybar;
                delete[] xbar;
                delete[] parent;
                delete[] cost;

                std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
                userTime += time_span.count();
            }
        }
    }
    catch (GRBException e)
    {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
    catch (...)
    {
        std::cout << "Error during callback" << std::endl;
    }
}