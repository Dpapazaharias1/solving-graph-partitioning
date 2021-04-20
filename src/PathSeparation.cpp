#include <algorithm>
#include <limits.h>
#include <random>
#include <chrono>
#include <queue>
#include "PathSeparation.h"

PathSeparation::PathSeparation(GRBEnv *env_, GRBVar **yvars, GRBVar *xvars, Graph &graph, int size, int p_)
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

                int rand_int = rand() % 100;

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
                int *parent = new int[n];
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
                    // -------- Tree Generation --------
                    
                    if(false && p > rand_int && G->total_weight == n) // Do not run for weighted instances
                    {
                        std::vector< std::vector<int> > treeAdj(n, std::vector<int>());
                        std::vector<int> numChildren(n, 0);
                        for(int i = 0; i < n; i++)
                        {
                            if (parent[i] > -1)
                            {
                                treeAdj[parent[i]].push_back(i);
                                numChildren[parent[i]]++;
                            }
                        }
                        treeCutLP(root, xbar, treeAdj, numChildren);
                    }

                    // -------- Branch number --------

                    if(false && p > rand_int)
                    {
                        std::vector< std::vector<int> > tree_adj(n, std::vector<int>());
                        std::vector<int> tree_degree(n, 0);
                        for(int i = 0; i < n; i++)
                        {
                            if( parent[i] > - 1)
                            {
                                tree_adj[i].push_back(parent[i]);
                                tree_adj[parent[i]].push_back(i);
                                tree_degree[i]++;
                                tree_degree[parent[i]]++;
                            }
                        }
                        compute_branch_number(xbar, tree_adj, tree_degree);
                    }
                    
                }

                

                // -------- Prim's Fractional Separation Subroutine ---------
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
                            treeWeight += xbar[j];
                            expr += x[j];
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
            int *parent = new int[n];
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

void PathSeparation::treeCutLP(int root, const double *xbar, const std::vector< std::vector<int> > &tree_adj, const std::vector<int> &num_children)
{
    GRBModel model = GRBModel(*env);

    std::vector<int> prefix_num_children(n, 0);
    std::partial_sum(std::begin(num_children), std::end(num_children), std::begin(prefix_num_children));
    
    // ---------- DFS Routine ------------
    int total_weight = 0;
    std::vector<bool> discovered(n, false);
    std::vector<int> subtree_weight(n, 0);
    std::vector<int> discover_weight(n, 0);
    compute_tree_weight(tree_adj, discovered, discover_weight, subtree_weight, total_weight, root);
    // ---------- DFS Routine ------------
    GRBVar zvar = model.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "z");
    
    //std::vector< std::vector<GRBVar> > hvar = (n, std::vector<GRBVar>());
    std::vector<GRBVar> hvar;
    std::vector< std::vector<GRBVar> > fvar(n, std::vector<GRBVar>());
    std::vector<  std::vector<GRBVar> > gvar(prefix_num_children.back(), std::vector<GRBVar>());
    std::vector<int> edge_costs;
    std::vector<int> edge_indices;
    int e = 0;
    for(int i = 0; i < n; i++)
    {
        for(int o = 0; o < G->r + 1; ++o)
        {
            fvar[i].push_back(model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
            if(o < G->Weight[i] || o > subtree_weight[i]) { fvar[i][o].set(GRB_DoubleAttr_LB, 2 * subtree_weight[root]); }  
        }
        if(num_children[i] == 0)
        {
            fvar[i][G->Weight[i]].set(GRB_DoubleAttr_UB, 0.0);       
        }
        for(int j = 0; j < num_children[i]; j++)
        {
            for(int k = G->EdgesBegin[i]; k < G->EdgesBegin[i] + G->Degree[i]; k++)
            {
                if(G->EdgeTo[k] == tree_adj[i][j])
                {
                    hvar.push_back(model.addVar(0.0, GRB_INFINITY, xbar[k], GRB_CONTINUOUS));
                    edge_costs.push_back(G->EdgeCost[k]);
                    edge_indices.push_back(k);
                    for(int o = 0; o < G->r + 1; ++o)
                    {
                        gvar[e].push_back(model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
                    }
                    e++;
                    break;
                }
            }   
        }
    }
    model.update();
    model.set(GRB_IntAttr_ModelSense, -1); // max
    model.set(GRB_IntParam_InfUnbdInfo, 1);
    model.set(GRB_IntParam_OutputFlag, 0);
    // Constraint 1: z >= f[root][o]
    for(int o = G->Weight[root]; o < G->r + 1; o++)
    {
        model.addConstr(zvar, GRB_LESS_EQUAL, fvar[root][o]);
    }

    // Constraint 2: f[i][o] = g[i][s][o] -> g[prefix[i - 1]+num_children[i]][o]

    for(int i = 0; i < n; i++)
    {
        if(num_children[i] > 0)
        {
            for(int o = G->Weight[i]; o < G->r; o++)
            {
                model.addConstr(fvar[i][o], GRB_EQUAL,  gvar[prefix_num_children[i] - 1][o]);
            }
            // Constraint 3: Cut
            for(int k = 1; k < G->r + 1; k++)
            {
                int first_idx = prefix_num_children[i] - num_children[i];
                int first_child = tree_adj[i].front();
                model.addConstr(gvar[first_idx][G->Weight[i]], GRB_LESS_EQUAL, edge_costs[first_idx] + fvar[first_child][k] + hvar[first_idx]);
            }
            // Constraint 4: Keep
            for(int o = G->Weight[i] + 1; o < G->r + 1; o++)
            {
                int first_idx = prefix_num_children[i] - num_children[i];
                int first_child = tree_adj[i].front();
                model.addConstr(gvar[first_idx][o], GRB_LESS_EQUAL, fvar[first_child][o - G->Weight[i]]);
            }
            
            for(int j = 1; j < num_children[i]; j++)
            {
                int edge_idx = prefix_num_children[i] - num_children[i] + j;
                int child = tree_adj[i][j];
                for(int o = 1; o < G->r + 1; o++)
                {
                    // Constraint 5: Cut
                    for(int k = 1; k < G->r + 1; k++)
                    {
                        model.addConstr(gvar[edge_idx][o], GRB_LESS_EQUAL, edge_costs[edge_idx] + gvar[edge_idx - 1][o] + fvar[child][k] + hvar[edge_idx]);
                    }
                    // Constraint 6: Keep
                    for(int k = 1; k < o + 1; k++)
                    {
                        model.addConstr(gvar[edge_idx][o], GRB_LESS_EQUAL, gvar[edge_idx - 1][k] + fvar[child][o - k]);
                    }
                }
            }
        }        
    }

    model.optimize();
    if(model.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
    {
        GRBLinExpr tree_cut = 0;
        double rhs = zvar.get(GRB_DoubleAttr_UnbdRay);
        for(int e = 0; e < n - 1; e++)
        {
            
            double coef = hvar[e].get(GRB_DoubleAttr_UnbdRay);
            if(coef > 0)
            {
                tree_cut += coef * x[edge_indices[e]];
            }
        }
        if(tree_cut.size() > 0 && rhs > 0)
        {
            lp_tree_cuts++;
            addCut(tree_cut, GRB_GREATER_EQUAL, rhs);
        }
    }
}

void PathSeparation::compute_branch_number(const double *xbar,
                                           std::vector< std::vector<int> > tree_adj, 
                                           std::vector<int> tree_degree)
{
    std::queue<int> Q;
    std::vector<int> tree_weight(n, 0);
    std::vector<int> branch_number(n, 0);
    for(int i = 0; i < n; i++){
        tree_weight[i] = G->Weight[i];
        if(tree_degree[i] == 1)
        {
            Q.push(i);
        }
    }
    int edge_param = std::min(G->r, G->total_weight - G->r);

    GRBLinExpr branching_cut = 0;
    double sum = 0;
    while(Q.size() > 1)
    {
        int u = Q.front();
        int v = tree_adj[u].front();
        for(int w : tree_adj[u])
        {
            if(tree_degree[w] > 0)
            {
                v = w;
                break;
            }
        }
        int b = std::min(tree_weight[u], G->total_weight - tree_weight[u]);
        int a = std::min(edge_param, b);
        for(int k = G->EdgesBegin[u]; k < G->EdgesBegin[u] + G->Degree[u]; k++)
        {
            if(G->EdgeTo[k] == v)
            {
                branching_cut += a * x[k];
                sum += a * xbar[k];
                break;
            }
        }

        
        
        tree_degree[u]--;
        tree_degree[v]--;
        tree_weight[v] += tree_weight[u];
        if(tree_degree[v] == 1)
        {
            Q.push(v);
        }
        Q.pop();
    }

    if(sum < G->total_weight - G->r)
    {
        branch_tree_cuts++;
        addCut(branching_cut, GRB_GREATER_EQUAL, G->total_weight - G->r);
    }
}

void PathSeparation::compute_tree_weight(const std::vector< std::vector<int> > &tree_adj,  
                                         std::vector<bool> &discovered, 
                                         std::vector<int> &s, 
                                         std::vector<int> &t,
                                         int &w, 
                                         int root)
{
    discovered[root] = true;
    s[root] = w;
    w += G->Weight[root];
    for(int v : tree_adj[root])
    {
        if(!discovered[v]) { compute_tree_weight(tree_adj, discovered, s, t, w, v); }
    }
    t[root] = w - s[root];
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