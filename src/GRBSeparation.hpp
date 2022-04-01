#ifndef GRBSEPARATION_H
#define GRBSEPARATION_H
#include "Graph.hpp"
#include "gurobi_c++.h"


class GPSeparation : public GRBCallback
{
    protected:
        GRBEnv env;
        std::vector<GRBVar> y;
        std::vector<GRBVar> x;
        Graph *G;
        int n;
        int m;
        int r;
        int p;
        int branch_tree_cuts;
        int lp_tree_cuts;
        float r_pct;
        const char *cut_type;
        int cut_code; // 0 - None, 1 - Prim, 2 - Knapsack, 3 - Tree Dynamic program
        std::vector<double> xsol;
        std::vector<double> ysol;

        void populate_x_MIPSOL()  {
            xsol = std::vector<double>(n * n, 0.0);
            for(int i = 0; i < n * n; i++) { xsol[i] = getSolution(x[i]); }
        }
        void populate_y_MIPSOL()  {
            ysol = std::vector<double>(2 * m, 0.0);
            for(int i = 0; i < 2 * m; i++) { ysol[i] = getSolution(y[i]); } 
        }
        void populate_x_MIPNODE() {
            xsol = std::vector<double>(n * n, 0.0);
            for(int i = 0; i < n * n; i++) { xsol[i] = getNodeRel(x[i]); }
        }
        void populate_y_MIPNODE() {
            ysol = std::vector<double>(2 * m, 0.0);
            for(int i = 0; i < 2 * m; i++) { ysol[i] = getNodeRel(y[i]); } 
        }

        int get_cut_code(const char* cut_type) {
            if(p == 0) {return 0;}
            if(strcmp(cut_type, "-prim") == 0) {return 1;}
            if(strcmp(cut_type, "-knap") == 0) {return 2;}
            if(strcmp(cut_type, "-tdp") == 0) {return 3;}
            return 0;
        };

        void connected_components(std::vector< std::vector<int>> &components, std::vector<int> &parent);
        void find_components(std::vector<bool> &processed, std::vector<int> &parent, std::vector<int> &component, int &component_weight, int root);
        void prim_heuristic(std::vector<double> &cost_to, std::vector<int> &parent, std::vector<double> &cost);
        void knapsack_cuts(std::vector<double> &cost_to, std::vector<int> &parent, std::vector<double> &cost);
        void tdp_separation(std::vector<double> &cost_to, std::vector<int> &parent, std::vector<double> &cost);
        void block_decomposition(std::vector<int> &tree_edges, const std::vector<int> &parent);
        void tdp_lp_separation(int root, const std::vector<double> &ybar, const std::vector< std::vector<int> > &tree_adj, const std::vector<int> &num_children, const std::vector<int> &tree_nodes);
        void compute_tree_weight(const std::vector<std::vector<int>> &tree_adj, std::vector<bool> &discovered, std::vector<int> &s, std::vector<int> &t, int &w, int root);
    public:
        int user_cuts;
        int lazy_cuts;
        double user_time;
        double lazy_time;

        GPSeparation(GRBEnv &env_, std::vector<GRBVar> &yvars, std::vector<GRBVar> &xvars, Graph &graph, int p_,int _r_pct, const char *_cut_type) : env(env_), y(yvars), x(xvars), G(&graph), p(p_), r_pct(_r_pct), cut_type(_cut_type)
        {
            n = graph.n;
            m = graph.m;
            r = graph.r;
            user_cuts = 0;
            lazy_cuts = 0;
            user_time = 0.0;
            lazy_time = 0.0;
            cut_code = get_cut_code(cut_type);
        }
};

class PathSeparation : public GPSeparation
{
    public:
        PathSeparation(
            GRBEnv env_, 
            std::vector<GRBVar> yvars, 
            std::vector<GRBVar> xvars, 
            Graph &graph, 
            int p_,
            float _r_pct, 
            const char *_cut_type
        ) : GPSeparation(env_, yvars, xvars, graph, p_, _r_pct, _cut_type) {};
    protected:
        void callback();
};

class TreeSeparation : public GPSeparation
{
    public:
        TreeSeparation(
            GRBEnv env_, 
            std::vector<GRBVar> yvars, 
            std::vector<GRBVar> xvars, 
            Graph &graph, 
            int p_, 
            float _r_pct, 
            const char *_cut_type
        ) : GPSeparation(env_, yvars, xvars, graph, p_, _r_pct, _cut_type) {};
    protected:
        void callback();
};


#endif