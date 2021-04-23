#include "Graph.h"
#include "PathSeparation.h"
#include "TreeSeparation.h"
#include "gurobi_c++.h"

void triangle(const Graph &G)
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
                if(i = j)
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
        std::cout << model.get(GRB_DoubleAttr_MIPGap) << " ";
        std::cout << model.get(GRB_DoubleAttr_Runtime) << " ";
        std::cout << model.get(GRB_DoubleAttr_NodeCount) << " ";
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

void flow(Graph &G, int p)
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
                if(i = j)
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
        model.set(GRB_IntParam_Cuts, 0);
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
            for(int j = 0; j < G.n; ++j)
            {
                constr += G.Weight[j] * x[(i * G.n) + j];
            }
            model.addConstr(constr, GRB_LESS_EQUAL, G.r - G.Weight[i]);
        }

        TreeSeparation cb = TreeSeparation(env, y, x, G, G.r, p);
        model.set(GRB_IntParam_PreCrush, 1);
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
        std::cout << model.get(GRB_DoubleAttr_NodeCount) << " ";
        
        // Branch and cut information cuts, callback time
        std::cout << p << " ";
        std::cout << cb.userCuts << " ";
        std::cout << cb.userTime << " ";
        // End line
        std::cout << std::endl;
    }
    catch (GRBException e)
    {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
}


void path(Graph &G, int p)
{

    // --------- Initialize Model and Environment ---------

    GRBEnv *env = new GRBEnv();
    GRBModel model = GRBModel(*env);

    GRBVar *x = new GRBVar[2 * G.m];
    GRBVar **y = new GRBVar *[G.n];

    int j;

    // --------- Initialize Variables ---------

    for (int i = 0; i < G.n; i++)
    {
        y[i] = new GRBVar[G.n];
    }

    std::string varName;

    // -------- x_e ---------

    for (int i = 0; i < G.n; i++)
    {
        for (int k = G.EdgesBegin[i]; k < G.EdgesBegin[i] + G.Degree[i]; k++)
        {
            j = G.EdgeTo[k];
            if (i < j)
            {
                varName = varName = "x_" + std::to_string(i) + "," + std::to_string(j);
                x[k] = model.addVar(0.0, 1.0, G.EdgeCost[k], GRB_BINARY, varName);
                for (int l = G.EdgesBegin[j]; l < G.EdgesBegin[j] + G.Degree[j]; l++)
                {
                    if (i == G.EdgeTo[l])
                    {
                        x[l] = x[k];
                    }
                }
            }
        }
    }

    // -------- y_i,j ---------

    for (int i = 0; i < G.n; i++)
    {
        for (int j = i + 1; j < G.n; j++)
        {
            varName = "y_" + std::to_string(i) + "," + std::to_string(j);
            y[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, varName);
            y[j][i] = y[i][j];
        }
        varName = "y_" + std::to_string(i) + "," + std::to_string(i);
        y[i][i] = model.addVar(0.0, 0.0, 0.0, GRB_CONTINUOUS, varName);
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
                model.addConstr(y[i][j] >= 1 - x[k]);
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
                expr += G.Weight[j] * y[i][j];
            }
        }
        model.addConstr(expr, GRB_LESS_EQUAL, G.r - G.Weight[i]);
    }

    // -------- Callback --------

    model.set(GRB_IntParam_PreCrush, 1);        // User Cuts
    model.set(GRB_IntParam_LazyConstraints, 1); // Lazy Cuts
    PathSeparation cb = PathSeparation(env, y, x, G, G.r, p);
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
    std::cout << cb.lazyCuts << " ";
    std::cout << cb.userCuts << " ";
    std::cout << cb.lazyTime << " ";
    std::cout << cb.userTime << " ";
    // Tree cut information
    std::cout << p << " ";
    std::cout << cb.branch_tree_cuts << " ";
    std::cout << cb.lp_tree_cuts << " ";
    // End line
    std::cout << std::endl;

    
    // -------- Memory Dellocation ---------

    for (int i = 0; i < G.n; i++)
    {
        delete[] y[i];
    }
    delete[] x;
    delete[] y;
    delete env;
}
