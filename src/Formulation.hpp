#ifndef _FORMULATION_H_
#define _FORMULATION_H_

#include "Graph.hpp"
#include "gurobi_c++.h"

void variable_setup(const Graph &G, GRBModel &model, std::vector<GRBVar> &x, std::vector<GRBVar> &y);
void triangle (const Graph &G, bool is_integer);
void flow(Graph &G, int p, float r_pct, const char* cut_type);
void path(Graph &G, int p, float r_pct, const char* cut_type);
void tree_cover_ip(Graph &G, const char* cut_type);
void tree_cover_formulation(Graph &G);
void dynamic_program_lp(Graph &G);
void tree_lp(Graph &G);
void flow_relax(Graph &G, int p, float rpct);
void compute_tree_weight(const std::vector< std::vector<int> > &tree_adj, const std::vector<int> &weight, std::vector<bool> &discovered,  std::vector<int> &s, std::vector<int> &t, int &w, int root);
#endif