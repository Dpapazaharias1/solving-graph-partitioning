#include <algorithm>
#include <limits.h>
#include <random>
#include <chrono>
#include "PathSeparation.h"

PathSeparation::PathSeparation(GRBVar **yvars, GRBVar *xvars, Graph &graph, int size)
{
    y = yvars;
    x = xvars;
    r = size;
    n = graph.n;
    m = graph.m;
    G = &graph;
    userCuts = 0;
    lazyCuts = 0;
    userTime = 0.0;
    lazyTime = 0.0;
}

void PathSeparation::callback()
{
    try
    {
        if (where == GRB_CB_MIPNODE)
        {
            if (getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
            {
                std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
                // -------- Retrieve Solution --------

                double **ybar = new double *[n];
                double *xbar = new double[2 * m];
                for (int i = 0; i < n; i++)
                {
                    ybar[i] = getNodeRel(y[i], n);
                }
                xbar = getNodeRel(x, 2 * m);

                // -------- Temp: Print Current Solution ---------

                //printSol(ybar, xbar);

                // -------- Place weight on edges ---------

                // We make the weight of using a cut-edge xbar_e

                std::vector<double> CostTo(2 * m, 1.0);
                for (int e = 0; e < 2 * m; e++)
                {
                    CostTo[e] = xbar[e];
                    // Some approximations of 0 are negative, causing negative cycles
                    if (xbar[e] < 0)
                    {
                        CostTo[e] = 0;
                    }
                }

                /* -------- Run Dijkstra's Fractional Separation Subroutine --------

                We will use Dijkstra's algorithm with each node 
                as the root to find violated path constraints.
                For each run of DA, we will keep a pointer to the array for parents and path cost
                The cut generation is as follows:
                    1. If a node's parent is -1 ignore, the node is not in the same component as root.
                    2. If a node's parent is root ignore, we have constraints for each edge.
                    3. If cost of path plus y_{ij} is less than 1 then we have a violated constraint, separate
                */

                int j;
                GRBLinExpr expr;
                double *parent = new double[n];
                double *cost = new double[n];

                for (int root = 0; root < n; root++)
                {
                    G->SP(CostTo, root, parent, cost);
                    for (int i = root + 1; i < n; i++) // We do not want to produce duplicate paths
                    {
                        /*
                            First condition: Violating Path
                            Second condition: We already have constraints for edges
                            Third condition: There is no path to nodes with parent -1
                        */

                        if (cost[i] + ybar[root][i] < 1 && parent[i] != root && parent[i] != -1)
                        {
                            expr = y[root][i]; // Constraint will be sum_Pij(x) + y_ij >= 1
                            j = i;
                            // Loop until we reach root (parent[root] = -1)
                            while (parent[j] != -1)
                            {
                                for (int e = G->EdgesBegin[j]; e < G->EdgesBegin[j] + G->Degree[j]; e++)
                                {
                                    if (parent[j] == G->EdgeTo[e])
                                    {
                                        expr += x[e]; // Need to get the edge # corresponding to i,j
                                        break;
                                    }
                                }
                                j = parent[j];
                            }
                            // Add cut, move to the next node
                            addCut(expr, GRB_GREATER_EQUAL, 1);
                            userCuts += 1;
                        }
                    }
                }

                // -------- Prim's Fractional Separation Subroutine ---------
                /*
                std::vector<int> treeNodes;
                srand((unsigned)time(0));
                int s = rand() % n;
                int treeWeight;
                treeNodes = G->Prim(CostTo, s, r, parent, cost);
                expr = 0;
                for (int i : treeNodes)
                {
                    if (parent[i] != -1)
                    {
                        j = parent[i];
                        treeWeight += xbar[j];
                        expr += x[j];
                    }
                }
                if (treeWeight < 1)
                {
                    addCut(expr, GRB_GREATER_EQUAL, 1);
                }
                */
                // -------- Memory Dellocation --------
                for (int i = 0; i < n; i++)
                {
                    delete[] ybar[i];
                }
                delete[] ybar;
                delete[] xbar;
                delete[] parent;
                delete[] cost;

                std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
                userTime += time_span.count();            
            }
        }

        if (where == GRB_CB_MIPSOL)
        {
            std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
            // -------- Retrieve Solution --------

            double **ybar = new double *[n];
            double *xbar = new double[2 * m];
            for (int i = 0; i < n; i++)
            {
                ybar[i] = getSolution(y[i], n);
            }
            xbar = getSolution(x, 2 * m);

            // -------- Temp: Print Current Solution ---------

            //printSol(ybar, xbar);

            // -------- Place weight on edges ---------

            // We make the weight of using a cut-edge xbar_e

            std::vector<double> CostTo(2 * m, 1.0);
            for (int e = 0; e < 2 * m; e++)
            {
                CostTo[e] = xbar[e];
                // Some approximations of 0 are negative, causing negative cycles
                if (xbar[e] < 0)
                {
                    CostTo[e] = 0;
                }
            }

            /* -------- Run Dijkstra's Integer Separation Subroutine --------

                We will use Dijkstra's algorithm with each node 
                as the root to find violated path constraints.
                For each run of DA, we will keep a pointer to the array for parents and path cost
                The cut generation is as follows:
                    1. If a node's parent is -1 ignore, the node is not in the same component as root.
                    2. If a node's parent is root ignore, we have constraints for each edge.
                    3. If the path has non-zero cost ignore, that path constraint is satisfied
                    4. Otherwise, use the parent list to get path to the root, add cut.
            */

            int j;
            GRBLinExpr expr;
            double *parent = new double[n];
            double *cost = new double[n];

            for (int root = 0; root < n; root++)
            {
                G->SP(CostTo, root, parent, cost);
                for (int i = root + 1; i < n; i++)
                {
                    /*
                        First condition: Violating Path
                        Second condition: We already have constraints for edges
                        Third condition: There is no path to nodes with parent -1
                    */

                    if (cost[i] == 0 && parent[i] != root && parent[i] != -1)
                    {
                        expr = y[root][i]; // Constraint will be sum_Pij(x) + y_ij >= 1
                        j = i;
                        // Loop until we reach root (parent[root] = -1)
                        while (parent[j] != -1)
                        {
                            for (int e = G->EdgesBegin[j]; e < G->EdgesBegin[j] + G->Degree[j]; e++)
                            {
                                if (parent[j] == G->EdgeTo[e])
                                {
                                    expr += x[e]; // Need to get the edge # corresponding to i,j
                                    break;
                                }
                            }
                            j = parent[j];
                        }
                        // Add cut, move to the next node
                        addLazy(expr, GRB_GREATER_EQUAL, 1);
                        lazyCuts += 1;
                    }
                }
            }
            // -------- Memory Dellocation --------

            for (int i = 0; i < n; i++)
            {
                delete[] ybar[i];
            }
            delete[] ybar;
            delete[] xbar;
            delete[] parent;
            delete[] cost;

            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
            lazyTime += time_span.count();            
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

void PathSeparation::printSol(double **ybar, double *xbar)
{
    std::cout << "----------- Connected Node Pairs -----------\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if (ybar[i][j] > 0.5)
            {
                std::cout << "(" << i << ", " << j << "): " << ybar[i][j] << "\n";
            }
        }
    }
    std::cout << "----------- Cut Edges -----------\n";
    int j;
    for (int i = 0; i < n; i++)
    {
        for (int k = G->EdgesBegin[i]; k < G->EdgesBegin[i] + G->Degree[i]; k++)
        {
            j = G->EdgeTo[k];
            if (i < j && xbar[k] > 0.5)
            {
                std::cout << "Cut edge: " << k << "\n";
            }
        }
    }
}