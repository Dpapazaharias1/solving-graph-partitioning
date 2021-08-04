#include "Graph.hpp"
#include "GRBSeparation.hpp"
#include "gurobi_c++.h"
#include <chrono>
#include <queue>
#include <numeric>
#include <cstring>

void triangle(const Graph &G, bool is_integer)
{
    try
    {
        GRBEnv *env = new GRBEnv();
        GRBModel model = GRBModel(*env);

        GRBVar *x = model.addVars(G.n * G.n, GRB_CONTINUOUS);
        GRBVar *y = model.addVars(2 * G.m, GRB_CONTINUOUS);
        std::cout << "xvars.." << std::endl;
        for(int i = 0; i < G.n; ++i)
        {
            for(int j = 0; j < G.n; ++j)
            {
                if(i == j){ x[(i * G.n) + j].set(GRB_DoubleAttr_UB, 0.0); }
                if(i < j) { x[(j * G.n) + i] = x[(i * G.n) + j]; }
            }
        }
        std::cout << "yvars.." << std::endl;
        for(int m = 0; m < 2 * G.m; ++m)
        {
            y[m].set(GRB_DoubleAttr_Obj, G.EdgeCost[m]);
            if(is_integer) {y[m].set(GRB_CharAttr_VType, GRB_BINARY);}
        }

        for (int i = 0; i < G.n; i++){
            for (int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++){
                int j = G.EdgeTo[k];
                if (i < j){
                    for (int l = G.EdgesBegin[j]; l < G.EdgesBegin[j] + G.Degree[j]; l++){
                        if (i == G.EdgeTo[l]){ y[l] = y[k]; }
                    }
                }
            }
        }

        model.set(GRB_IntAttr_ModelSense, 1);
        model.set(GRB_IntParam_OutputFlag, 0);
        model.set(GRB_DoubleParam_TimeLimit, 7200);
        model.update();
        std::cout << "constr1.." << std::endl;
        for(int i = 0; i < G.n; ++i)
        {
            for(int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; ++k)
            {
                int j = G.EdgeTo[k];
                if(i < j)
                {
                    if (is_integer) {
                        x[(i * G.n) + j].set(GRB_CharAttr_VType, GRB_BINARY);
                        x[(j * G.n) + i].set(GRB_CharAttr_VType, GRB_BINARY);
                    }
                    model.addConstr(y[k], GRB_EQUAL, 1 - x[(i * G.n) + j]);
                    model.addConstr(y[k], GRB_EQUAL, 1 - x[(j * G.n) + i]);
                }
            }
        }
        std::cout << "constr2.." << std::endl;
        for(int i = 0; i < G.n; ++i)
        {
            for(int j = i + 1; j < G.n; ++j)
            {
                for(int k = j + 1; k < G.n; ++k)
                {
                    model.addConstr(x[(i * G.n) + j] + x[(i * G.n) + k] - x[(j * G.n) + k], GRB_LESS_EQUAL, 1.0);
                    model.addConstr(x[(i * G.n) + j] - x[(i * G.n) + k] + x[(j * G.n) + k], GRB_LESS_EQUAL, 1.0);
                    model.addConstr(-x[(i * G.n) + j] + x[(i * G.n) + k] + x[(j * G.n) + k], GRB_LESS_EQUAL, 1.0);
                }
            }
        }
        GRBLinExpr constr;
        for(int i = 0; i < G.n; ++i)
        {
            constr = 0;
            for(int j = 0; j < G.n; ++j)
            {
                constr += G.Weight[j] * x[(i * G.n) + j];
            }
            model.addConstr(constr, GRB_LESS_EQUAL, G.r - G.Weight[i]);
        }

        model.optimize();
        // General solver information graphtype n m r best_obj gap% runtime nodecount
        std::cout << G.graphtype << " ";
        std::cout << G.n << " ";
        std::cout << G.m << " ";
        std::cout << G.r << " ";
        std::cout << model.get(GRB_DoubleAttr_ObjVal) << " ";
        if (is_integer) { std::cout << model.get(GRB_DoubleAttr_MIPGap) << " "; }
        std::cout << model.get(GRB_DoubleAttr_Runtime) << " ";
        if (is_integer) { std::cout << model.get(GRB_DoubleAttr_NodeCount) << " "; }
        
        // Branch and cut information cuts, callback time - NOT USED HERE

        // End line
        std::cout << std::endl;
    }
    catch (GRBException e)
    {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
}

void flow(Graph &G, int p, float r_pct, const char* cut_type)
{
    try
    {
        GRBEnv *env = new GRBEnv();
        GRBModel model = GRBModel(*env);

        GRBVar *x = model.addVars(G.n * G.n, GRB_CONTINUOUS);
        GRBVar *y = model.addVars(2 * G.m, GRB_BINARY);
        std::cout << "xvars.." << std::endl;
        for(int i = 0; i < G.n; ++i)
        {
            for(int j = 0; j < G.n; ++j)
            {
                if(i == j)
                {
                    x[(i * G.n) + j].set(GRB_DoubleAttr_UB, 0.0);
                }
                if(i < j)
                {
                    x[(j * G.n) + i] = x[(i * G.n) + j];
                }
            }
        }
        std::cout << "yvars.." << std::endl;
        for(int m = 0; m < 2 * G.m; ++m)
        {
            y[m].set(GRB_DoubleAttr_Obj, G.EdgeCost[m]); // TODO: Make it work for general costs~[1, 1000]
            
        }

        for (int i = 0; i < G.n; i++){
            for (int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++){
                int j = G.EdgeTo[k];
                if (i < j){
                    for (int l = G.EdgesBegin[j]; l < G.EdgesBegin[j] + G.Degree[j]; l++){
                        if (i == G.EdgeTo[l]){ y[l] = y[k]; }
                    }
                }
            }
        }

        model.set(GRB_IntAttr_ModelSense, 1);
        model.set(GRB_IntParam_OutputFlag, 0);
        model.set(GRB_DoubleParam_TimeLimit, 7200);
        model.update();
        std::cout << "constr1.." << std::endl;
        for(int i = 0; i < G.n; ++i)
        {
            for(int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; ++k)
            {
                int j = G.EdgeTo[k];
                if(i < j)
                {
                    x[(i * G.n) + j].set(GRB_CharAttr_VType, GRB_BINARY);
                    x[(j * G.n) + i].set(GRB_CharAttr_VType, GRB_BINARY);
                    model.addConstr(y[k], GRB_EQUAL, 1 - x[(i * G.n) + j]);
                    model.addConstr(y[k], GRB_EQUAL, 1 - x[(j * G.n) + i]);
                }
            }
        }
        std::cout << "constr2.." << std::endl;
        for(int i = 0; i < G.n; ++i){
            for(int q = G.EdgesBegin[i]; q < G.EdgesBegin[i] + G.Degree[i]; q++){
                int j = G.EdgeTo[q];
                for(int k = 0; k < G.n; ++k){
                    if (i < j && k != i && k != j)
                    {
                        model.addConstr(x[(i * G.n) + j] + x[(i * G.n) + k] - x[(j * G.n) + k], GRB_LESS_EQUAL, 1.0);
                        model.addConstr(x[(i * G.n) + j] - x[(i * G.n) + k] + x[(j * G.n) + k], GRB_LESS_EQUAL, 1.0);
                        model.addConstr(-x[(i * G.n) + j] + x[(i * G.n) + k] + x[(j * G.n) + k], GRB_LESS_EQUAL, 1.0);
                    }
                }
            }
        }

        GRBLinExpr constr;
        for(int i = 0; i < G.n; ++i)
        {
            constr = 0;
            for(int j = 0; j < G.n; ++j) { constr += G.Weight[j] * x[(i * G.n) + j]; }
            model.addConstr(constr, GRB_LESS_EQUAL, G.r - G.Weight[i]);
        }

        TreeSeparation cb = TreeSeparation(env, y, x, G, p, r_pct, cut_type);
        if(p > 0 && (strcmp(cut_type, "-prim") == 0 || strcmp(cut_type, "-knap") == 0 || strcmp(cut_type, "-tdp") == 0))
        {
            model.set(GRB_IntParam_Cuts, 0);
            model.set(GRB_IntParam_PreCrush, 1);
            model.setCallback(&cb);
        }

        model.optimize();
        // General solver information graphtype n m r best_obj gap% runtime nodecount
        std::cout << G.graphtype << " ";
        std::cout << G.n << " ";
        std::cout << G.m << " ";
        std::cout << G.r << " ";
        std::cout << model.get(GRB_DoubleAttr_ObjVal) << " ";
        std::cout << model.get(GRB_DoubleAttr_MIPGap) << " ";
        std::cout << model.get(GRB_DoubleAttr_Runtime) << " ";
        std::cout << model.get(GRB_DoubleAttr_NodeCount) << " ";
        
        // Branch and cut information cuts, callback time
        std::cout << cut_type << " ";
        std::cout << p << " ";
        std::cout << r_pct << " ";
        std::cout << cb.user_cuts << " ";
        std::cout << cb.user_time << " ";
        // End line
        std::cout << std::endl;
    }
    catch (GRBException e)
    {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
}


void path(Graph &G, int p, float r_pct, const char* cut_type)
{
    GRBEnv *env = new GRBEnv();
    GRBModel model = GRBModel(*env);

    GRBVar *y = new GRBVar[2 * G.m];
    GRBVar *x = new GRBVar[G.n * G.n];

    int j;

    std::string varName;

    for (int i = 0; i < G.n; i++) {
        for (int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++) {
            j = G.EdgeTo[k];
            if (i < j) {
                varName = varName = "y_" + std::to_string(i) + "," + std::to_string(j);
                y[k] = model.addVar(0.0, 1.0, G.EdgeCost[k], GRB_BINARY, varName);
                for (int l = G.EdgesBegin[j]; l < G.EdgesBegin[j] + G.Degree[j]; l++) {
                    if (i == G.EdgeTo[l]) { y[l] = y[k]; }
                }
            }
        }
    }

    // -------- y_i,j ---------

    for (int i = 0; i < G.n; i++)
    {
        for (int j = i + 1; j < G.n; j++)
        {
            varName = "x_" + std::to_string(i) + "," + std::to_string(j);
            x[(i * G.n) + j] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, varName);
            x[(j * G.n) + i] = x[(i * G.n) + j];
        }
        varName = "x_" + std::to_string(i) + "," + std::to_string(i);
        x[(i * G.n) + i] = model.addVar(0.0, 0.0, 0.0, GRB_CONTINUOUS, varName);
    }
    model.set(GRB_DoubleParam_TimeLimit, 7200);
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_IntAttr_ModelSense, 1);
    model.update();

    // -------- Constraints ---------

    for (int i = 0; i < G.n; i++)
    {
        for (int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++)
        {
            j = G.EdgeTo[k];
            if (i < j)
            {
                model.addConstr(x[(i * G.n) + j]>= 1 - y[k]);
            }
        }
    }

    GRBLinExpr expr = 0;
    for (int i = 0; i < G.n; i++)
    {
        expr = 0;
        for (int j = 0; j < G.n; j++)
        {
            if (i != j)
            {
                expr += G.Weight[j] * x[(i * G.n) + j];
            }
        }
        model.addConstr(expr, GRB_LESS_EQUAL, G.r - G.Weight[i]);
    }

    // -------- Callback --------

    model.set(GRB_IntParam_PreCrush, 1);        // User Cuts
    model.set(GRB_IntParam_LazyConstraints, 1); // Lazy Cuts
    PathSeparation cb = PathSeparation(env, y, x, G, p, r_pct, cut_type);
    if(p > 0) {model.set(GRB_IntParam_Cuts, 0);}
    model.setCallback(&cb);

    // -------- Optimize Model --------

    model.optimize();

    // -------- Prepare Output --------

    // General solver information graphtype n m r best_obj gap% runtime nodecount
    std::cout << G.graphtype << " ";
    std::cout << G.n << " ";
    std::cout << G.m << " ";
    std::cout << G.r << " ";
    std::cout << model.get(GRB_DoubleAttr_ObjVal) << " ";
    std::cout << model.get(GRB_DoubleAttr_MIPGap) << " ";
    std::cout << model.get(GRB_DoubleAttr_Runtime) << " ";
    std::cout << model.get(GRB_DoubleAttr_NodeCount) << " ";
    // Branch and cut information cuts, callback time
    std::cout << cb.lazy_cuts << " ";
    std::cout << cb.user_cuts << " ";
    std::cout << cb.lazy_time << " ";
    std::cout << cb.user_time << " ";
    // Tree cut information
    std::cout << cut_type << " ";
    std::cout << p << " ";
    std::cout << r_pct << " ";
    //std::cout << cb.branch_tree_cuts << " ";
    //std::cout << cb.lp_tree_cuts << " ";
    // End line
    std::cout << std::endl;

    
    // -------- Memory Dellocation ---------

    delete[] x;
    delete[] y;
    delete env;
}

void tree_cover_ip(Graph &G, const char* cut_type)
{
    // --------- Initialize Model and Environment ---------
        std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
        GRBEnv *env = new GRBEnv();
        GRBModel model = GRBModel(*env);

        GRBVar *y = new GRBVar[2 * G.m];
        GRBVar *x = new GRBVar[G.n * G.n];
        int userCuts = 0;
        

        for (int i = 0; i < G.n; i++){
            for (int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++){
                int j = G.EdgeTo[k];
                if (i < j){
                    y[k] = model.addVar(0.0, 1.0, G.EdgeCost[k], GRB_BINARY);
                    for (int l = G.EdgesBegin[j]; l < G.EdgesBegin[j] + G.Degree[j]; l++){
                        if (i == G.EdgeTo[l]){ y[l] = y[k]; }
                    }
                }
            }
        }
        
        for(int i = 0; i < G.n * G.n; i++) {x[i] = model.addVar(0.0, 0.0, 0.0, GRB_CONTINUOUS);}

        model.set(GRB_DoubleParam_TimeLimit, 7200);
        model.set(GRB_IntAttr_ModelSense, 1);
        model.set(GRB_IntParam_OutputFlag, 0);
        model.set(GRB_IntParam_LazyConstraints, 1);
        model.update();
        
        TreeSeparation cb = TreeSeparation(env, y, x, G, 100, 1.00, cut_type);
        
        model.setCallback(&cb);
        model.optimize();

        // General solver information graphtype n m r best_obj gap% runtime nodecount
        std::cout << G.graphtype << " ";
        std::cout << G.n << " ";
        std::cout << G.m << " ";
        std::cout << G.r << " ";
        std::cout << model.get(GRB_DoubleAttr_ObjVal) << " ";
        std::cout << model.get(GRB_DoubleAttr_MIPGap) << " ";
        std::cout << model.get(GRB_DoubleAttr_Runtime) << " ";
        // Branch and cut information cuts, callback time
        std::cout << cut_type << " ";
        std::cout << cb.lazy_cuts << " ";
        std::cout << cb.lazy_time << " ";
        // Tree cut information
        // End line
        std::cout << std::endl;
}


void tree_cover_formulation(Graph &G)
{
    try {
        // --------- Initialize Model and Environment ---------
        std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
        GRBEnv *env = new GRBEnv();
        GRBModel model = GRBModel(*env);

        GRBVar *y = new GRBVar[2 * G.m];
        int userCuts = 0;
        

        for (int i = 0; i < G.n; i++){
            for (int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++){
                int j = G.EdgeTo[k];
                if (i < j){
                    y[k] = model.addVar(0.0, 1.0, G.EdgeCost[k], GRB_CONTINUOUS);
                    for (int l = G.EdgesBegin[j]; l < G.EdgesBegin[j] + G.Degree[j]; l++){
                        if (i == G.EdgeTo[l]){ y[l] = y[k]; }
                    }
                }
            }
        }

        model.set(GRB_DoubleParam_TimeLimit, 7200);
        model.set(GRB_IntAttr_ModelSense, 1);
        model.set(GRB_IntParam_OutputFlag, 1);
        model.set(GRB_IntParam_OutputFlag, 0);
        //model.set(GRB_IntParam_PreCrush, 1);
        model.update();
        /*
        TreeSeparation cb = TreeSeparation(env, y, x, G, G.r, 0);
        cb.is_tcf = true;
        model.setCallback(&cb);
        */

        model.optimize();

        GRBModel MST_model = GRBModel(*env);

        GRBVar *z = new GRBVar[G.n * (G.n + 1)];
        GRBVar *f = new GRBVar[G.n * (G.n + 1)];

        for (int i = 0; i < G.n; i++){
            for(int j = 0; j < G.n + 1; j++)
            {
                std::string z_name = "z_" + std::to_string(i) + "," + std::to_string(j);
                std::string f_name = "f_" + std::to_string(i) + "," + std::to_string(j);
                z[(i * (G.n + 1)) + j] = MST_model.addVar(0.0, 0.0, 0.0, GRB_BINARY, z_name);
                f[(i * (G.n + 1)) + j] = MST_model.addVar(0.0, 0.0, 0.0, GRB_CONTINUOUS, f_name);
            }
        }
        for (int i = 0; i < G.n; i++)
        {
            for(int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++)
            {
                int j = G.EdgeTo[k];
                if(i < j)
                {
                    z[(i * (G.n + 1)) + j].set(GRB_DoubleAttr_Obj, y[k].get(GRB_DoubleAttr_X));
                    z[(i * (G.n + 1)) + j].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                }
                f[(i * (G.n + 1)) + j].set(GRB_DoubleAttr_UB, GRB_INFINITY);
                
            }
            z[(i * (G.n + 1)) + G.n].set(GRB_DoubleAttr_UB, GRB_INFINITY);
            f[(i * (G.n + 1)) + G.n].set(GRB_DoubleAttr_UB, 1.0);
        }

        MST_model.set(GRB_IntAttr_ModelSense, 1);
        MST_model.set(GRB_IntParam_OutputFlag, 0);
        MST_model.update();

        for (int i = 0; i < G.n; i++){
            for(int j = 0; j < G.n + 1; j++){
                MST_model.addConstr(f[(i * (G.n + 1)) + j], GRB_LESS_EQUAL, z[(i * (G.n + 1)) + j] * G.r);
            }
        }


        GRBConstr* flow_constraint = new GRBConstr[G.n];
        GRBLinExpr expr;
        for(int i = 0; i < G.n; i++){
            expr = 0;
            expr += f[(i * (G.n + 1)) + G.n];
            for(int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++)
            {
                int j = G.EdgeTo[k];
                expr += f[(i * (G.n + 1)) + j];
                if(j < G.n) { expr -= f[(j * (G.n + 1)) + i]; }
            }
            flow_constraint[i] = MST_model.addConstr(expr, GRB_EQUAL, 0.0);
        }
        expr = 0;
        for(int i = 0; i < G.n; i++)
        {
            expr += f[(i * (G.n + 1)) + G.n];
        }
        MST_model.addConstr(expr, GRB_EQUAL, G.r);
        
        expr = 0;
        for(int i = 0; i < G.n; i++) {
            for(int j = 0; j < G.n; j++) { expr += z[(i * (G.n + 1)) + j]; }
        }
        MST_model.addConstr(expr, GRB_EQUAL, G.r);

        GRBLinExpr tree_cut;
        bool violated = true;
        while(violated)
        {
            violated = false;
            model.optimize();
            for (int i = 0; i < G.n; i++){
                for(int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++){
                    int j = G.EdgeTo[k];
                    z[(i * (G.n + 1)) + j].set(GRB_DoubleAttr_Obj, y[k].get(GRB_DoubleAttr_X));
                }
            }

            for(int i = 0; i < G.n; i++)
            {
                f[(i * (G.n + 1)) + G.n].set(GRB_DoubleAttr_UB, G.r);
                z[(i * (G.n + 1)) + G.n].set(GRB_DoubleAttr_Obj, 2.0 * G.n);
                flow_constraint[i].set(GRB_DoubleAttr_RHS, G.r);

                MST_model.update();
                MST_model.optimize();
                tree_cut = 0;
                if(MST_model.get(GRB_DoubleAttr_ObjVal) < 1 - 1e-5) {
                    violated = true;
                    //std::cout << "Tree at root " << i << " objective: " << MST_model.get(GRB_DoubleAttr_ObjVal) << std::endl;
                    //std::cout << "Edges: ";
                    for(int j = 0; j < G.n; j++) {
                        for(int k = G.EdgesBegin[j]; k < G.EdgesBegin[j] + G.Degree[j]; k++) {
                            int l = G.EdgeTo[k];
                            if(z[(j * (G.n + 1)) + l].get(GRB_DoubleAttr_X) > 0.5)
                            {
                                tree_cut += y[k];
                            }
                        }
                    }
                    //std::cout << std::endl;
                    model.addConstr(tree_cut, GRB_GREATER_EQUAL, 1.0);
                    userCuts++;
                }

                f[(i * (G.n + 1)) + G.n].set(GRB_DoubleAttr_UB, 1);
                z[(i * (G.n + 1)) + G.n].set(GRB_DoubleAttr_Obj, 0);
                flow_constraint[i].set(GRB_DoubleAttr_RHS, 0);
            }
        }
        
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> userTime = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
        // General solver information graphtype n m r best_obj gap% runtime nodecount
        std::cout << G.graphtype << " ";
        std::cout << G.n << " ";
        std::cout << G.m << " ";
        std::cout << G.r << " ";
        std::cout << model.get(GRB_DoubleAttr_ObjVal) << " ";
        std::cout << model.get(GRB_DoubleAttr_Runtime) << " ";
        // Branch and cut information cuts, callback time
        std::cout << userCuts << " ";
        std::cout << userTime.count() << " ";
        // Tree cut information
        // End line
        std::cout << std::endl;
    }
    catch (GRBException e)
    {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
}

void compute_tree_weight(const std::vector< std::vector<int> > &tree_adj,
                         const std::vector<int> &weight,
                        std::vector<bool> &discovered, 
                        std::vector<int> &s, 
                        std::vector<int> &t,
                        int &w, 
                        int root)
{
    discovered[root] = true;
    s[root] = w;
    w += weight[root];
    for(int v : tree_adj[root])
    {
        if(!discovered[v]) { compute_tree_weight(tree_adj, weight, discovered, s, t, w, v); }
    }
    t[root] = w - s[root];
}




void dynamic_program_lp(Graph &G)
{
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    GRBEnv *env = new GRBEnv();
    GRBModel model = GRBModel(*env);

    GRBVar *y = new GRBVar[2 * G.m];
    int userCuts = 0;
    for (int i = 0; i < G.n; i++){
        for (int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++){
            int j = G.EdgeTo[k];
            if (i < j){
                y[k] = model.addVar(0.0, 1.0, G.EdgeCost[k], GRB_CONTINUOUS);
                for (int l = G.EdgesBegin[j]; l < G.EdgesBegin[j] + G.Degree[j]; l++){
                    if (i == G.EdgeTo[l]){ y[l] = y[k]; }
                }
            }
        }
    }

    model.set(GRB_DoubleParam_TimeLimit, 7200);
    model.set(GRB_IntAttr_ModelSense, 1);
    model.set(GRB_IntParam_OutputFlag, 1);
    model.set(GRB_IntParam_OutputFlag, 0);
    //model.set(GRB_IntParam_PreCrush, 1);
    model.update();

    std::vector<std::vector<int>> tree_adj(G.n, std::vector<int>());
    std::vector<int> num_children(G.n, 0);
    std::vector<int> parent(G.n, -1);
    std::vector<bool> discovered(G.n, false);
    std::queue<int> Q;
    discovered[0] = true;
    Q.push(0);
    while(!Q.empty())
    {
        int i = Q.front();
        Q.pop();
        
        for(int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++)
        {
            int j = G.EdgeTo[k];
            if(!discovered[j])
            {
                discovered[j] = true;
                Q.push(j);
                tree_adj[i].push_back(j);
                num_children[i]++;
                parent[j] = i;
            }
        }
    }

    std::vector<int> prefix_num_children(G.n, 0);
    std::partial_sum(std::begin(num_children), std::end(num_children), std::begin(prefix_num_children));

    for(int i = 0; i < G.n; i++){ discovered[i] = false; }
    std::vector<int> subtree_weight(G.n, 0);
    std::vector<int> discover_weight(G.n, 0);
    int total_weight = 0;
    compute_tree_weight(tree_adj, G.Weight, discovered, discover_weight, subtree_weight, total_weight, 0);



    GRBModel tree_dp = GRBModel(*env);
    

    GRBVar zvar = tree_dp.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "z");
    std::vector<GRBVar> hvar;
    std::vector< std::vector<GRBVar> > fvar(G.n, std::vector<GRBVar>());
    std::vector<  std::vector<GRBVar> > gvar(prefix_num_children.back(), std::vector<GRBVar>());
    std::vector<int> edge_indices;
    int e = 0;
    for(int i = 0; i < G.n; i++)
    {
        for(int o = 0; o < G.r + 1; ++o)
        {
            fvar[i].push_back(tree_dp.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
            if(o < G.Weight[i] || o > subtree_weight[i]) { fvar[i][o].set(GRB_DoubleAttr_LB, 2 * subtree_weight[0]); }
        }
        if(num_children[i] == 0) { fvar[i][G.Weight[i]].set(GRB_DoubleAttr_UB, 0.0); }
        for(int j = 0; j < num_children[i]; j++)
        {
            for(int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++)
            {
                if(G.EdgeTo[k] == tree_adj[i][j])
                {
                    hvar.push_back(tree_dp.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
                    edge_indices.push_back(k);
                    for(int o = 0; o < G.r + 1; ++o) {gvar[e].push_back(tree_dp.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));}
                    e++;
                    break;
                }
            }
        }
    }

    tree_dp.update();
    tree_dp.set(GRB_IntAttr_ModelSense, -1);
    tree_dp.set(GRB_IntParam_InfUnbdInfo, 1);
    tree_dp.set(GRB_IntParam_OutputFlag, 0);

    for(int o = G.Weight[0]; o < G.r + 1; o++) { tree_dp.addConstr(zvar, GRB_LESS_EQUAL, fvar[0][o]); }

    for(int i = 0; i < G.n; i++)
    {
        if(num_children[i] > 0)
        {
            int first_idx = prefix_num_children[i] - num_children[i];
            int first_child = tree_adj[i].front();
            for(int o = G.Weight[i]; o < G.r + 1; o++) {tree_dp.addConstr(fvar[i][o], GRB_EQUAL, gvar[prefix_num_children[i] - 1][o]);}
            for(int k = 1; k < G.r + 1; k++) {tree_dp.addConstr(gvar[first_idx][G.Weight[i]], GRB_LESS_EQUAL, G.EdgeCost[edge_indices[first_idx]] + fvar[first_child][k] + hvar[first_idx]);}
            for(int o = G.Weight[i] + 1; o < G.r + 1; o++) {tree_dp.addConstr(gvar[first_idx][o], GRB_LESS_EQUAL, fvar[first_child][o - G.Weight[i]]);}

            for(int j = 1; j < num_children[i]; j++)
            {
                int edge_idx = first_idx + j;
                int child = tree_adj[i][j];
                for(int o = 1; o < G.r + 1; o++)
                {
                    for(int k = 1; k < G.r + 1; k++)
                    {
                        tree_dp.addConstr(gvar[edge_idx][o], GRB_LESS_EQUAL, G.EdgeCost[edge_indices[edge_idx]]+ gvar[edge_idx - 1][o] + fvar[child][k] + hvar[edge_idx]);
                    }
                    for(int k = 1; k < o + 1; k++)
                    {
                        tree_dp.addConstr(gvar[edge_idx][o], GRB_LESS_EQUAL, gvar[edge_idx - 1][k] + fvar[child][o - k]);
                    }
                }
            }
        }
    }

    model.optimize();
    bool violated = true;
    GRBLinExpr tree_cut;
    while (violated) {
        violated = false;
        for(int k = 0; k < G.n - 1; k++) { hvar[k].set(GRB_DoubleAttr_Obj, -y[edge_indices[k]].get(GRB_DoubleAttr_X));}
        tree_dp.update();
        tree_dp.optimize();
        if(tree_dp.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        {
            std::cout << "Optimal solution: " << tree_dp.get(GRB_DoubleAttr_ObjVal) << std::endl;
        }
        if(tree_dp.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
        {
            std::cout << "Model is infeasible" << std::endl;
        }

        if(tree_dp.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
        {
            violated = true;
            tree_cut = 0;
            double rhs = zvar.get(GRB_DoubleAttr_UnbdRay);
            
            for(int k = 0; k < G.n - 1; k++)
            {
                double coef = hvar[k].get(GRB_DoubleAttr_UnbdRay);
                if(coef > 0) { tree_cut += coef * y[edge_indices[k]];}    
            }
            model.addConstr(tree_cut, GRB_GREATER_EQUAL, rhs);
            userCuts++;
        }
        model.optimize();
    }
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> userTime = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
    // General solver information graphtype n m r best_obj gap% runtime nodecount
    std::cout << G.graphtype << " ";
    std::cout << G.n << " ";
    std::cout << G.m << " ";
    std::cout << G.r << " ";
    std::cout << model.get(GRB_DoubleAttr_ObjVal) << " ";
    std::cout << model.get(GRB_DoubleAttr_Runtime) << " ";
    // Branch and cut information cuts, callback time
    std::cout << userCuts << " ";
    std::cout << userTime.count() << " ";
    // Tree cut information
    // End line
    std::cout << std::endl;
}

void tree_lp(Graph &G)
{
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    try{
        GRBEnv *env = new GRBEnv();
        GRBModel model = GRBModel(*env);

        GRBVar *y = new GRBVar[2 * G.m];
        int userCuts = 0;
        std::cout << "y variables" << std::endl;
        for (int i = 0; i < G.n; i++){
            for (int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++){
                int j = G.EdgeTo[k];
                if (i < j){
                    y[k] = model.addVar(0.0, 1.0, G.EdgeCost[k], GRB_CONTINUOUS);
                    for (int l = G.EdgesBegin[j]; l < G.EdgesBegin[j] + G.Degree[j]; l++){
                        if (i == G.EdgeTo[l]){ y[l] = y[k]; }
                    }
                }
            }
        }

        model.set(GRB_DoubleParam_TimeLimit, 7200);
        model.set(GRB_IntAttr_ModelSense, 1);
        model.set(GRB_IntParam_OutputFlag, 1);
        model.set(GRB_IntParam_OutputFlag, 0);
        model.update();

        std::cout << "getting rooted out tree" << std::endl;

        std::vector<std::vector<int>> tree_adj(G.n, std::vector<int>());
        std::vector<int> num_children(G.n, 0);
        std::vector<int> parent(G.n, -1);
        std::vector<bool> discovered(G.n, false);
        std::queue<int> Q;
        discovered[0] = true;
        Q.push(0);
        while(!Q.empty())
        {
            int i = Q.front();
            Q.pop();
            
            for(int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++)
            {
                int j = G.EdgeTo[k];
                if(!discovered[j])
                {
                    discovered[j] = true;
                    Q.push(j);
                    tree_adj[i].push_back(j);
                    num_children[i]++;
                    parent[j] = i;
                }
            }
        }
        
        std::vector<int> prefix_num_children(G.n, 0);
        std::partial_sum(std::begin(num_children), std::end(num_children), std::begin(prefix_num_children));

        for(int i = 0; i < G.n; i++){ discovered[i] = false; }
        std::vector<int> subtree_weight(G.n, 0);
        std::vector<int> discover_weight(G.n, 0);
        int total_weight = 0;
        compute_tree_weight(tree_adj, G.Weight, discovered, discover_weight, subtree_weight, total_weight, 0);

        std::vector<GRBVar> alpha(G.r + 1);
        std::vector<std::vector<GRBVar> > beta(G.n, std::vector<GRBVar>());
        std::vector< std::vector< std::vector<GRBVar> > > gamma(G.n - 1);
        std::vector< std::vector< std::vector<GRBVar> > > rho(G.n - 1);
        std::vector<int> edge_number;
        std::string vname;
        std::cout << "Creating variables: alpha" << std::endl;
        for(int o = 1; o < G.r + 1; o++)
        {
            vname = "alpha_" + std::to_string(o);
            alpha[o] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, vname);
        }
        std::cout << "Creating variables: beta" << std::endl;
        for(int i = 0; i < G.n; i++){
            for(int o = 0; o < G.r + 1; o++)
            {
                vname = "beta_" + std::to_string(i) + "," + std::to_string(o);
                beta[i].push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, vname));
            }
        }
        std::cout << "Creating variables: gamma & rho" << std::endl;
        int e = 0;
        for(int i = 0; i < G.n; i++)
        {
            for(int j = 0; j < num_children[i]; j++)
            {
                for(int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++)
                {
                    if(G.EdgeTo[k] == tree_adj[i][j]){ edge_number.push_back(k);}
                }
                gamma[e] = std::vector<std::vector<GRBVar> >(G.r + 1);
                rho[e] = std::vector<std::vector<GRBVar> >(G.r + 1);
                
                for(int o = 0; o < G.r + 1; o++){
                    gamma[e][o] = std::vector<GRBVar>();
                    rho[e][o] = std::vector<GRBVar>();
                    for(int k = 0; k < G.r + 1; k++)
                    {
                        vname = "gamma_" + std::to_string(e) + "," + std::to_string(o) + "," + std::to_string(k);
                        gamma[e][o].push_back(model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, vname));
                    }
                    for(int k = 0; k < o; k++)
                    {
                        vname = "rho_" + std::to_string(e) + "," + std::to_string(o) + "," + std::to_string(k);
                        rho[e][o].push_back(model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, vname));
                    }
                }
                e++;
            }
        }
        

        model.update();
        model.set(GRB_IntAttr_ModelSense, 1);

        std::cout << "Edge cut" << std::endl;

        GRBLinExpr edge_cut;
        for(int edge = 0; edge < G.n - 1; edge++){
            edge_cut = 0;
            for(int o = 1; o < G.r + 1; o++) {
                for(int k = 1; k < G.r + 1; k++){edge_cut += gamma[edge][o][k];}
            }
            model.addConstr(edge_cut, GRB_LESS_EQUAL, y[edge_number[edge]]);
        }
        std::cout << "alpha cut" << std::endl;
        GRBLinExpr expr = 0;
        for(int o = 1; o < G.r + 1; o++)
        {
            expr += alpha[o];
        }
        model.addConstr(expr, GRB_GREATER_EQUAL, 1.0);
        std::cout << "beta cut" << std::endl;
        for(int o = 1; o < G.r + 1; o++)
        {
            model.addConstr(-alpha[o] + beta[0][o], GRB_GREATER_EQUAL, 0.0);
        }
        
        std::cout << "gamma cuts" << std::endl;
        GRBLinExpr expr1,expr2, expr3;
        for(int i = 0; i < G.n; i++) {
            for(int j = 0; j < num_children[i]; j++) {
                e = prefix_num_children[i] - num_children[i] + j;
                for(int q = 1; q < G.r + 1; q++) {
                    expr1 = 0;
                    expr1 += beta[tree_adj[i][j]][q];
                    for(int o = 1; o < G.r + 1; o++)
                    {
                        expr1 += -gamma[e][o][q]; 
                    }
                    for(int o = 1; o < G.r - q + 1; o++){ 
                        expr1 += -rho[e][o + q][o];
                    }
                    
                    model.addConstr(expr1, GRB_GREATER_EQUAL, 0.0);
                }
                
                if(j == num_children[i] - 1) {
                    for(int o = 1; o < G.r + 1; o++){
                        expr2 = 0;
                        expr2 += -beta[i][o];
                        for(int q = 1; q < G.r + 1; q++) { expr2 += gamma[e][o][q];}
                        for(int q = 1; q < o; q++) { expr2 += rho[e][o][o - q]; }
                        model.addConstr(expr2, GRB_GREATER_EQUAL, 0.0);
                    }
                }
                if(j > 0  && e > 0) {
                    for(int o = 1; o < G.r + 1; o++) {
                        expr3 = 0;
                        for(int q = 1; q < G.r + 1; q++) {
                            expr3 += gamma[e - 1][o][q];
                            expr3 += -gamma[e][o][q];
                            if (q < o){
                                expr3 += rho[e - 1][o][o - q];
                            }
                            if (q < G.r + 1 - o) {
                                expr3 += -rho[e][o + q][o];
                            }
                        }
                        model.addConstr(expr3, GRB_GREATER_EQUAL, 0.0);
                    }
                }
            }
        }
        

        model.optimize();
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> userTime = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
        // General solver information graphtype n m r best_obj gap% runtime nodecount
        std::cout << G.graphtype << " ";
        std::cout << G.n << " ";
        std::cout << G.m << " ";
        std::cout << G.r << " ";
        std::cout << model.get(GRB_DoubleAttr_ObjVal) << " ";
        std::cout << model.get(GRB_DoubleAttr_Runtime) << " ";
        std::cout << std::endl;
    }
    catch (GRBException e)
    {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
}


void flow_relax(Graph &G, int p, float rpct)
{
    try
    {
        std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
        GRBEnv *env = new GRBEnv();
        GRBModel model = GRBModel(*env);

        GRBVar *x = model.addVars(G.n * G.n, GRB_CONTINUOUS);
        GRBVar *y = model.addVars(2 * G.m, GRB_CONTINUOUS);
        std::cout << "xvars.." << std::endl;
        for(int i = 0; i < G.n; ++i)
        {
            for(int j = 0; j < G.n; ++j)
            {
                if(i == j)
                {
                    x[(i * G.n) + j].set(GRB_DoubleAttr_UB, 0.0);
                }
                if(i < j)
                {
                    x[(j * G.n) + i] = x[(i * G.n) + j];
                }
            }
        }
        std::cout << "yvars.." << std::endl;
        for(int m = 0; m < 2 * G.m; ++m)
        {
            y[m].set(GRB_DoubleAttr_UB, 1.0);
            y[m].set(GRB_DoubleAttr_Obj, G.EdgeCost[m]); // TODO: Make it work for general costs~[1, 1000]
            
        }

        for (int i = 0; i < G.n; i++){
            for (int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++){
                int j = G.EdgeTo[k];
                if (i < j){
                    for (int l = G.EdgesBegin[j]; l < G.EdgesBegin[j] + G.Degree[j]; l++){
                        if (i == G.EdgeTo[l]){ y[l] = y[k]; }
                    }
                }
            }
        }

        model.set(GRB_IntAttr_ModelSense, 1);
        model.set(GRB_IntParam_OutputFlag, 0);
        model.set(GRB_DoubleParam_TimeLimit, 7200);
        model.update();
        std::cout << "constr1.." << std::endl;
        for(int i = 0; i < G.n; ++i)
        {
            for(int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; ++k)
            {
                int j = G.EdgeTo[k];
                if(i < j)
                {
                    x[(i * G.n) + j].set(GRB_DoubleAttr_UB, 1.0);
                    x[(j * G.n) + i].set(GRB_DoubleAttr_UB, 1.0);
                    x[(i * G.n) + j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                    x[(j * G.n) + i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                    model.addConstr(y[k], GRB_EQUAL, 1 - x[(i * G.n) + j]);
                    model.addConstr(y[k], GRB_EQUAL, 1 - x[(j * G.n) + i]);
                }
            }
        }
        std::cout << "constr2.." << std::endl;
        for(int i = 0; i < G.n; ++i){
            for(int q = G.EdgesBegin[i]; q < G.EdgesBegin[i] + G.Degree[i]; q++){
                int j = G.EdgeTo[q];
                for(int k = 0; k < G.n; ++k){
                    if (i < j && k != i && k != j)
                    {
                        model.addConstr(x[(i * G.n) + j] + x[(i * G.n) + k] - x[(j * G.n) + k], GRB_LESS_EQUAL, 1.0);
                        model.addConstr(x[(i * G.n) + j] - x[(i * G.n) + k] + x[(j * G.n) + k], GRB_LESS_EQUAL, 1.0);
                        model.addConstr(-x[(i * G.n) + j] + x[(i * G.n) + k] + x[(j * G.n) + k], GRB_LESS_EQUAL, 1.0);
                    }
                }
            }
        }

        GRBLinExpr constr;
        for(int i = 0; i < G.n; ++i)
        {
            constr = 0;
            for(int j = 0; j < G.n; ++j)
            {
                constr += G.Weight[j] * x[(i * G.n) + j];
            }
            model.addConstr(constr, GRB_LESS_EQUAL, G.r - G.Weight[i]);
        }

        
        // While LOOP
        bool violated_cut = true;
        std::vector<int> tree_nodes;
        std::vector<double> cost_to(2 * G.m, 0.0);
        int user_cuts = 0;
        int root;
        double tree_weight;
        int node_weight;
        int tree_limit = G.r * rpct;
        int *parent = new int[G.n];
        double *cost = new double[G.n];
        std::vector< std::vector<int> > tree_adj(G.n, std::vector<int>());
        std::vector<int> tree_degree(G.n, 0);
        std::vector<int> prefix_children(G.n, 0);
        double current_time = 0;
	while(violated_cut && p > 0 && current_time <= 7200.00)
        {
            violated_cut = false;
            GRBModel tree_dp = GRBModel(*env);
            std::vector<GRBVar> hvar;
            std::vector< std::vector<GRBVar> > fvar;
            std::vector<  std::vector<GRBVar> > gvar;
            model.optimize();
            std::cout << "Current LP Value: "  << model.get(GRB_DoubleAttr_ObjVal) << std::endl;
            for(root = 0; root < G.n; root++)
            {
                //root = rand() % G.n;
                node_weight = 0;
                tree_weight = 0.0;
                tree_nodes.clear();
                for(int m = 0; m < 2 * G.m; m++) {cost_to[m] = y[m].get(GRB_DoubleAttr_X);}
                tree_nodes = G.Prim(cost_to, root, tree_limit, parent, cost);
                std::fill(tree_degree.begin(), tree_degree.end(), 0);
                std::fill(tree_adj.begin(), tree_adj.end(), std::vector<int>());
                for(int i : tree_nodes)
                {   if (parent[i] > -1)
                    {
                        tree_adj[G.EdgeTo[parent[i]]].push_back(i);
                        tree_degree[G.EdgeTo[parent[i]]]++;
                    }
                }

                int total_weight = 0;
                std::vector<bool> discovered(G.n, false);
                std::vector<int> subtree_weight(G.n, 0);
                std::vector<int> discover_weight(G.n, 0);
                compute_tree_weight(tree_adj, G.Weight, discovered, discover_weight, subtree_weight, total_weight, root);

                //std::partial_sum(std::begin(tree_degree), std::end(tree_degree), std::begin(prefix_children));
                fvar = std::vector< std::vector<GRBVar> >(G.n, std::vector<GRBVar>());
                gvar = std::vector< std::vector<GRBVar> >(tree_nodes.size() - 1, std::vector<GRBVar>());
                std::vector<int> edge_begin(G.n, -1);
                std::vector<int> edge_costs;
                std::vector<int> edge_index;
                int e = 0;
                GRBVar zvar = tree_dp.addVar(0.0, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
                for(int i : tree_nodes)
                {
                    for(int o = 0; o < G.r + 1; ++o)
                    {
                        std::string fname = "f_" + std::to_string(i) + ","  + std::to_string(o);
                        
                        fvar[i].push_back(tree_dp.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, fname));
                        
                        if(o < G.Weight[i] || o > subtree_weight[i]) { fvar[i][o].set(GRB_DoubleAttr_LB, 2 * subtree_weight[root]); }
                        
                    }
                    
                    if(tree_degree[i] == 0) {fvar[i][G.Weight[i]].set(GRB_DoubleAttr_UB, 0.0);}
                    if(tree_degree[i] > 0 ) { edge_begin[i] = e; }
                    
                    for(int j = 0; j < tree_degree[i]; j++)
                    {   
                        for(int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++)
                        {
                            if(G.EdgeTo[k] == tree_adj[i][j]){
                                hvar.push_back(tree_dp.addVar(0.0, GRB_INFINITY, -cost_to[k], GRB_CONTINUOUS));
                                edge_costs.push_back(cost_to[k]);
                                edge_index.push_back(k);
                                for(int o = 0; o < G.r + 1; o++)
                                {
                                    gvar[e].push_back(tree_dp.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS));
                                }
                                e++;
                                break;
                            }
                        }
                    }
                }
                tree_dp.update();
                tree_dp.set(GRB_IntAttr_ModelSense, -1);
                tree_dp.set(GRB_IntParam_InfUnbdInfo, 1);
                tree_dp.set(GRB_IntParam_OutputFlag, 0);
                for(int o = G.Weight[root]; o < G.r + 1; o++) {tree_dp.addConstr(zvar, GRB_LESS_EQUAL, fvar[root][o]);}

                
                for(int i : tree_nodes) {
                    if(tree_degree[i] > 0) {
                        int first_idx = edge_begin[i];
                        int first_child = tree_adj[i].front();
                        for(int o = G.Weight[i]; o < G.r + 1; o++) { tree_dp.addConstr(fvar[i][o], GRB_EQUAL,  gvar[first_idx + tree_degree[i] - 1][o]);}
                        
                        for(int k = 1; k < G.r + 1; k++)
                        {
                            tree_dp.addConstr(gvar[first_idx][G.Weight[i]], GRB_LESS_EQUAL, edge_costs[first_idx] + fvar[first_child][k] + hvar[first_idx]);
                        }
                        
                        for(int o = G.Weight[i] + 1; o < G.r + 1; o++)
                        {
                            tree_dp.addConstr(gvar[first_idx][o], GRB_LESS_EQUAL, fvar[first_child][o - G.Weight[i]]);
                        }
                        for(int j = 1; j < tree_degree[i]; j++)
                        {
                            int edge_idx = edge_begin[i] + j;
                            int child = tree_adj[i][j];
                            for(int o = 1; o < G.r + 1; o++)
                            {
                                for(int k = 1; k < G.r + 1; k++)
                                {
                                    tree_dp.addConstr(gvar[edge_idx][o], GRB_LESS_EQUAL, edge_costs[edge_idx] + gvar[edge_idx - 1][o] + fvar[child][k] + hvar[edge_idx]);
                                }
                                for(int k = 1; k < o + 1; k++)
                                {
                                    tree_dp.addConstr(gvar[edge_idx][o], GRB_LESS_EQUAL, gvar[edge_idx - 1][k] + fvar[child][o - k]);
                                }
                            }
                        }
                    }
                }
                while(true)
                {
                    tree_dp.optimize();
                    //model.write("treedp.lp");
                    if(tree_dp.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
                    {
                        //std::cout << "No cut" << std::endl;
                        break;
                    }
                    if(tree_dp.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
                    {
                        //std::cout << "Model Infeasible" << std::endl;
                        break;
                    }
                    if(tree_dp.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
                    {
                        GRBLinExpr tree_cut = 0;
                        double rhs = zvar.get(GRB_DoubleAttr_UnbdRay);
                        for(int e = 0; e < tree_nodes.size() - 1; e++)
                        {
                            
                            double coef = hvar[e].get(GRB_DoubleAttr_UnbdRay);
                            if(coef > 0)
                            {
                                tree_cut += coef * y[edge_index[e]];
                            }
                        }
                        if(tree_cut.size() > 0 && rhs > 0)
                        {
                            violated_cut = true;
                            user_cuts++;
                            std::cout << "new cut" << tree_nodes.size() << std::endl;
                            model.addConstr(tree_cut, GRB_GREATER_EQUAL, rhs);
                        }
                        
                    }
                    //break;
                    model.optimize();
                    for(int e = 0; e < tree_nodes.size() - 1; e++)
                    {
                        hvar[e].set(GRB_DoubleAttr_Obj, -y[edge_index[e]].get(GRB_DoubleAttr_X));
                    }
                    std::cout << "(Old tree) Current LP Value: "  << model.get(GRB_DoubleAttr_ObjVal) ", cuts added: " << user_cuts << std::endl;
                }
                fvar.clear();
                gvar.clear();
                hvar.clear();
        	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        	std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
		current_time = time_span.count();
            } // end while
        }
        if (p == 0) {model.optimize();}
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast< std::chrono::duration<double> >(t1 - t0);
        // General solver information graphtype n m r best_obj gap% runtime nodecount
        std::cout << G.graphtype << " ";
        std::cout << G.n << " ";
        std::cout << G.m << " ";
        std::cout << G.r << " ";
        std::cout << model.get(GRB_DoubleAttr_ObjVal) << " ";
        std::cout << time_span.count() << " ";
        
        // Branch and cut information cuts, callback time
        std::cout << p << " ";
        std::cout << user_cuts << " ";
       // std::cout << cb.userTime << " ";
        // End line
        std::cout << std::endl;
    }
    catch (GRBException e)
    {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
}
