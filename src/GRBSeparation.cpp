#include "GRBSeparation.hpp"
#include <chrono>
#include <queue>
#include <cstring>
#include <numeric>
#include <limits.h>

void GPSeparation::prim_heuristic(std::vector<double> &cost_to, std::vector<int> &parent, std::vector<double> &cost)
{
    GRBLinExpr expr;
    std::vector<int> tree_nodes;
    srand((unsigned) time(0));
    int s = rand() % n;
    double tree_weight = 0.0;
    int node_weight = 0;
    G->Prim(cost_to, s, r, parent, cost, tree_nodes);
    int j;
    expr = 0;
    for(int i : tree_nodes)
    {
        node_weight += G->Weight[i];
        if(parent[i] != -1)
        {
            j = parent[i];
            tree_weight += ysol[j];
            expr += y[j];
        }
    }
    if (tree_weight < 1.0 && node_weight > r)
    {
        addCut(expr, GRB_GREATER_EQUAL, 1.0);
        user_cuts++;
    }
}


void GPSeparation::knapsack_cuts(std::vector<double> &cost_to, std::vector<int> &parent, std::vector<double> &cost)
{
    GRBLinExpr expr;
    std::vector<int> tree_nodes;
    srand((unsigned) time(0));
    int s = rand() % n;
    double tree_weight = 0.0;
    int node_weight = 0;
    int tree_limit = r * r_pct;
    G->Prim(cost_to, s, tree_limit, parent, cost, tree_nodes);

    std::vector< std::vector<int> > tree_adj(n, std::vector<int>());
    std::vector<int> tree_degree(n, 0);
    for(int i : tree_nodes)
    {
        node_weight += G->Weight[i];
        if(parent[i] != -1)
        {
            tree_adj[i].push_back(G->EdgeTo[parent[i]]);
            tree_adj[G->EdgeTo[parent[i]]].push_back(i);
            tree_degree[i]++;
            tree_degree[G->EdgeTo[parent[i]]]++;
        }
    }

    std::queue<int> Q;
    std::vector<int> subtree_weight(n, 0);
    std::vector<int> branch_number(n, 0);
    int total_weight = 0;
    for(int i = 0; i < n; i++)
    {
        subtree_weight[i] = G->Weight[i];
        if (tree_degree[i] > 0) {total_weight += G->Weight[i];}
        if (tree_degree[i] == 1) {Q.push(i);}
    }
    int edge_param = std::min(r, total_weight - r);
    double sum = 0;
    expr = 0;
    while (Q.size() > 1)
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
        int b = std::min(subtree_weight[u], total_weight - subtree_weight[u]);
        int a = std::min(edge_param, b);
        for(int k = G->EdgesBegin[u]; k < G->EdgesBegin[u] + G->Degree[u]; k++)
        {
            if(G->EdgeTo[k] == v)
            {
                expr += a * y[k];
                sum += a * ysol[k];
                break;
            }
        }
        tree_degree[u]--;
        tree_degree[v]--;
        subtree_weight[v] += subtree_weight[u];
        if(tree_degree[v] == 1) {Q.push(v);}
        Q.pop();
    }
    if (sum < total_weight - r)
    {
        user_cuts++;
        addCut(expr, GRB_GREATER_EQUAL, total_weight - r);
    }
}

void GPSeparation::tdp_separation(std::vector<double> &cost_to, std::vector<int> &parent, std::vector<double> &cost)
{
    std::vector<int> tree_nodes;
    srand((unsigned) time(0));
    int s = rand() % n;
    double tree_weight = 0.0;
    int node_weight = 0;
    int tree_limit = r * r_pct;
    G->Prim(cost_to, s, tree_limit, parent, cost, tree_nodes);

    std::vector< std::vector<int> > tree_adj(n, std::vector<int>());
    std::vector<int> tree_degree(n, 0);
    for(int i : tree_nodes)
    {
        if(parent[i] > - 1){
            tree_adj[G->EdgeTo[parent[i]]].push_back(i);
            tree_degree[G->EdgeTo[parent[i]]]++;
        }
    }

    tdp_lp_separation(s, ysol, tree_adj, tree_degree, tree_nodes);
}

void GPSeparation::tdp_lp_separation(int root, 
                               const std::vector<double> &ybar, 
                               const std::vector< std::vector<int> > &tree_adj, 
                               const std::vector<int> &num_children, 
                               const std::vector<int> &tree_nodes)
{
    GRBModel model = GRBModel(env);

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
    for(int i : tree_nodes)
    {
        for(int o = 0; o < G->r + 1; ++o)
        {
            std::string fname = "f_" + std::to_string(i) + ","  + std::to_string(o);
            fvar[i].push_back(model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, fname));
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
                    std::string hname = "h_" + std::to_string(i) + ","  + std::to_string(tree_adj[i][j]);
                    hvar.push_back(model.addVar(0.0, GRB_INFINITY, -ybar[k], GRB_CONTINUOUS, hname));
                    edge_costs.push_back(G->EdgeCost[k]);
                    edge_indices.push_back(k);
                    for(int o = 0; o < G->r + 1; ++o)
                    {
                        std::string gname = "g_" + std::to_string(i) + ","  + std::to_string(tree_adj[i][j]) +"," + std::to_string(o);
                        gvar[e].push_back(model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, gname));
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

    for(int i : tree_nodes)
    {
        if(num_children[i] > 0)
        {
            for(int o = G->Weight[i]; o < G->r + 1; o++)
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
        std::cout << "No cut" << std::endl;
    }
    if(model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
    {
        std::cout << "Model Infeasible" << std::endl;
    }
    if(model.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
    {
        GRBLinExpr tree_cut = 0;
        double rhs = zvar.get(GRB_DoubleAttr_UnbdRay);
        for(int e = 0; e < tree_nodes.size() - 1; e++)
        {
            
            double coef = hvar[e].get(GRB_DoubleAttr_UnbdRay);
            if(coef > 0)
            {
                tree_cut += coef * y[edge_indices[e]];
            }
        }
        if(tree_cut.size() > 0 && rhs > 0)
        {
            user_cuts++;
            std::cout << "new cut" << std::endl;
            addCut(tree_cut, GRB_GREATER_EQUAL, rhs);
        }
    }
}

void GPSeparation::compute_tree_weight(const std::vector< std::vector<int> > &tree_adj,  
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


void GPSeparation::block_decomposition(std::vector<int> &tree_edges, const std::vector<int> &parent) 
{
    std::vector<float> edge_coef(2 * m, 0.0);
    std::vector< std::vector<int> > tree_adj(n, std::vector<int>());
    std::vector< std::vector<int> > tree_edge_list(n, std::vector<int>());
    std::vector<int> tree_degree(n, 0);
    int root;
    //std::cout << "Tree Edges: ";
    for(int e : tree_edges)
    {
        edge_coef[e] = 1.0;
        int j = G->EdgeTo[e];
        int i = parent[j];
        //std::cout << "(" << i << ", " << j << "), ";
        if (parent[i] < 0) {root = i;}
        tree_edge_list[i].push_back(e);
        tree_edge_list[j].push_back(e);
        tree_adj[i].push_back(j);
        tree_adj[j].push_back(i);
        tree_degree[i]++;
        tree_degree[j]++;
    }
    //std::cout << std::endl;
    //std::cout << "Tree nodes: ";
    std::vector<int> tree_nodes;
    for(int i = 0; i < n; i++) {
        if(tree_degree[i] > 0) {
            //std::cout << i << " ";
            tree_nodes.push_back(i);
        }
    }
    //std::cout << std::endl;

    std::vector<int> processed(n, false);
    std::vector<int> block;
    bool found_cycle;
    for(int i : tree_nodes)
    {
        if(!processed[i])
        {
            block.clear();
            found_cycle = false;
            for(int j : tree_adj[i])
            {
                if(!processed[j]) {block.push_back(j);}
            }
            //check if neighbors have an edge
            for(int j = 0 ; j < block.size(); j++)
            {
                for(int k = j + 1; k < block.size(); k++)
                {
                    if (found_cycle) {break;}
                    int u = block[j];
                    int v = block[k];
                    for(int e = G->EdgesBegin[u]; e < G->EdgesBegin[u] + G->Degree[u]; e++)
                    {
                        if(G->EdgeTo[e] == v) {
                            //std::cout << "1. Add edge (" << u << ", " << v  << "), creating cycle with node " << i << std::endl;
                            found_cycle = true;
                            tree_edges.push_back(e);
                            edge_coef[e] = 0.5;
                            // Now find coefficients of edge (i, u) and (i, v)
                            for(int s : tree_edge_list[i])
                            {
                                if (G->EdgeTo[s] == u)
                                {
                                    edge_coef[s] = 0.5;
                                }
                                if (G->EdgeTo[s] == v)
                                {
                                    edge_coef[s] = 0.5;
                                }
                                if (G->EdgeTo[s] == i)
                                {
                                    if (s >= G->EdgesBegin[u] && s < G->EdgesBegin[u] + G->Degree[u]) {edge_coef[s] = 0.5;}
                                    if (s >= G->EdgesBegin[v] && s < G->EdgesBegin[v] + G->Degree[v]) {edge_coef[s] = 0.5;}
                                }
                            }
                            break;
                        }
                    }
                }
                processed[j] = true;
            }
            if (!found_cycle) {
                for(int j : block)
                {
                    for (int k : tree_adj[j])
                    {
                        if(found_cycle) { processed[k] = true; }
                        else{
                            if(k == i || processed[k]) {continue;}
                            processed[k] = true;
                            for(int e = G->EdgesBegin[k]; e < G->EdgesBegin[k] + G->Degree[k]; e++)
                            {
                                if(G->EdgeTo[e] == i) {
                                    //std::cout << "2. Add edge (" << i << ", " << k  << "), creating cycle with node " << j << std::endl;
                                    found_cycle = true;
                                    tree_edges.push_back(e);
                                    edge_coef[e] = 0.5;
                                    // Now find coefficients of edge (i, j) and (j, k)
                                    for(int s : tree_edge_list[i])
                                    {
                                        if (G->EdgeTo[s] == j)
                                        {
                                            edge_coef[s] = 0.5;
                                        }
                                        if (G->EdgeTo[s] == i)
                                        {
                                            if (s >= G->EdgesBegin[j] && s < G->EdgesBegin[j] + G->Degree[j]) {edge_coef[s] = 0.5;}
                                        }
                                    }
                                    for(int s : tree_edge_list[j])
                                    {
                                        if (G->EdgeTo[s] == k)
                                        {
                                            edge_coef[s] = 0.5;
                                        }
                                        if (G->EdgeTo[s] == j)
                                        {
                                            if (s >= G->EdgesBegin[k] && s < G->EdgesBegin[k] + G->Degree[k]) {edge_coef[s] = 0.5;}
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            else {
                for(int j : block)
                {
                    for(int k : tree_adj[j]) {processed[k] = true;}
                }
            }
        }
    }

    GRBLinExpr expr = 0;
    for(int e : tree_edges)
    {
        expr += edge_coef[e] * y[e];
    }
    addLazy(expr, GRB_GREATER_EQUAL, 1.0);
    
}

void GPSeparation::connected_components(std::vector< std::vector<int>> &components, std::vector<int> &parent)
{   
    std::vector<int> component;
    int component_weight = 0;
    std::vector<bool> processed(n, false);
    for(int i = 0; i < n; i++)
    {
        if(!processed[i])
        {
            component_weight = 0;
            component.clear();
            find_components(processed, parent, component, component_weight, i);
            if (component_weight > r)
            {
                components.push_back(component);
            }
        }
    }
}

void GPSeparation::find_components(
    std::vector<bool> &processed,
    std::vector<int> &parent,
    std::vector<int> &component, 
    int &component_weight, 
    int root
)
{
    processed[root] = true;
    component_weight += G->Weight[root];
    for(int e = G->EdgesBegin[root]; e < G->EdgesBegin[root] + G->Degree[root]; e++)
    {
        if (!processed[G->EdgeTo[e]] && ysol[e] < 0.5 && component_weight < r + 1) // if child has not been processed and edge is not in cut
        {
            component.push_back(e);
            parent[G->EdgeTo[e]] = root;
            find_components(processed, parent, component, component_weight, G->EdgeTo[e]);
        }
    }
    
}

void PathSeparation::callback()
{
    try
    {
        if(where == GRB_CB_MIPNODE)
        {
            if(getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
            {
                std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
                populate_x_MIPNODE();
                populate_y_MIPNODE();

                std::vector<double> cost_to(2 * m, 1.0);
                for(int e = 0; e < 2 * m; e++)
                {
                    cost_to[e] = ysol[e];
                    if(ysol[e] < 0) { cost_to[e] = 0.0; }
                }

                int j;
                GRBLinExpr expr;
                std::vector<int> parent(n, -1);
                std::vector<double> cost(n, INT_MAX);

                for (int root = 0; root < n; root++) {
                    G->SP(cost_to, root, parent, cost);
                    for(int i = root + 1; i < n; i++) {
                        if (cost[i] + xsol[(root * n) + i] < 0.999 && parent[i] != root && parent[i] != -1) 
                        {
                            expr = x[(root * n) + i];
                            j = i;
                            while (parent[j] != -1) {
                                for(int e = G->EdgesBegin[j]; e < G->EdgesBegin[j] + G->Degree[j]; e++) {
                                    if(parent[j] == G->EdgeTo[e]) {
                                        expr += y[e];
                                        break;
                                    }
                                }
                                j = parent[j];
                            } // while - path reconstruction 
                            addCut(expr, GRB_GREATER_EQUAL, 1.0);
                            user_cuts++;
                        } // if cut is violated
                    }
                } // Shortest Path Separation

                // If running PATH+
                srand(time(NULL));
                int rand_int = rand() % 100;

                if(cut_code == 1 && p > rand_int) {
                    prim_heuristic(cost_to, parent, cost);
                }
                if(cut_code == 2 && p > rand_int) {
                    knapsack_cuts(cost_to, parent, cost);
                }
                if(cut_code == 3 && p > rand_int) {
                    tdp_separation(cost_to, parent, cost);
                }

                std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
                user_time += time_span.count();
                
            } // if MIPNODE - Optimal
        } // if MIPNODE

        if(where == GRB_CB_MIPSOL)
        {
            std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
            populate_x_MIPSOL();
            populate_y_MIPSOL();

            std::vector<double> cost_to(2 * m, 1.0);
            for(int e = 0; e < 2 * m; e++)
            {
                cost_to[e] = ysol[e];
                if(ysol[e] < 0) { cost_to[e] = 0.0; }
            }

            int j;
            GRBLinExpr expr;
            std::vector<int> parent(n, -1);
            std::vector<double> cost(n, INT_MAX);
            

            for (int root = 0; root < n; root++) {
                G->SP(cost_to, root, parent, cost);
                for(int i = root + 1; i < n; i++) {
                    if (cost[i] + xsol[(root * n) + i] < 0.999 && parent[i] != root && parent[i] != -1) 
                    {
                        expr = x[(root * n) + i];
                        j = i;
                        while (parent[j] != -1) {
                            for(int e = G->EdgesBegin[j]; e < G->EdgesBegin[j] + G->Degree[j]; e++) {
                                if(parent[j] == G->EdgeTo[e]) {
                                    expr += y[e];
                                    break;
                                }
                            }
                            j = parent[j];
                        } // while - path reconstruction
                        addLazy(expr, GRB_GREATER_EQUAL, 1.0);
                        lazy_cuts++;
                    } // if cut is violated
                }
            } // Shortest Path Separation

            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
            lazy_time += time_span.count();
        } // if MIPSOL
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

void TreeSeparation::callback()
{
    try
    {
        if(where == GRB_CB_MIPNODE)
        {
            if(getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
            {
                std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
                
                srand(time(NULL));
                int rand_int = rand() % 100;

                populate_x_MIPNODE();
                populate_y_MIPNODE();

                std::vector<double> cost_to(2 * m, 1.0);
                for(int e = 0; e < 2 * m; e++)
                {
                    cost_to[e] = ysol[e];
                    if(ysol[e] < 0) { cost_to[e] = 0.0; }
                }

                std::vector<int> parent(n, -1);
                std::vector<double> cost(n, INT_MAX);

                if(cut_code == 1 && p > rand_int) {
                    prim_heuristic(cost_to, parent, cost);
                }
                if(cut_code == 2 && p > rand_int) {
                    knapsack_cuts(cost_to, parent, cost);
                }
                if(cut_code == 3 && p > rand_int) {
                    tdp_separation(cost_to, parent, cost);
                }

                std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
                user_time += time_span.count();
                
            } // if MIPNODE - Optimal
        } // if MIPNODE

        if(where == GRB_CB_MIPSOL)
        {
            std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
            populate_x_MIPSOL();
            populate_y_MIPSOL();
            //std::cout << "Starting cut" << std::endl;
            std::vector< std::vector<int> > cover_trees;
            std::vector<int> parent(n, -1);
            connected_components(cover_trees, parent);

            if(strcmp(cut_type, "-tcf") == 0) {
                GRBLinExpr expr;
                //std::cout << "Violating components " << cover_trees.size() << std::endl; 
                for(std::vector<int> tree_edges : cover_trees)
                {
                    expr = 0;
                    for (int e : tree_edges) { 
                        expr += y[e];
                    }
                    addLazy(expr, GRB_GREATER_EQUAL, 1.0);
                    lazy_cuts++;
                }
            }

            if(strcmp(cut_type, "-block") == 0) {
                GRBLinExpr expr; 
                for(std::vector<int> tree_edges : cover_trees)
                {
                    block_decomposition(tree_edges, parent);
                    lazy_cuts++;
                }
            }

            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
            lazy_time += time_span.count();
        } // if MIPSOL
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