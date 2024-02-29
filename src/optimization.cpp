/************************************************************************
*
* ADOBE CONFIDENTIAL
* ___________________
*
* Copyright [first year code created] Adobe
* All Rights Reserved.
*
* NOTICE: All information contained herein is, and remains
* the property of Adobe and its suppliers, if any. The intellectual
* and technical concepts contained herein are proprietary to Adobe
* and its suppliers and are protected by all applicable intellectual
* property laws, including trade secret and copyright laws.
* Dissemination of this information or reproduction of this material
* is strictly forbidden unless prior written permission is obtained
* from Adobe.
*************************************************************************
*/

#include "optimization.h"

Eigen::VectorXd solve_QP_with_ineq(Eigen::SparseMatrix<double> Q, Eigen::VectorXd g, Eigen::VectorXd x_0, 
                                   Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b){
                                        // Gurobi calling code borrowed from CoMISo

    //----------------------------------------------
    // 0. set up gurobi environment
    //----------------------------------------------
    GRBEnv env     = GRBEnv();
    GRBModel model = GRBModel(env);

    double time_limit = 1e10;
    model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);

    //----------------------------------------------
    // 1. allocate variables
    //----------------------------------------------
    // GUROBI variables

    size_t n_vars = x_0.size(),
           n_fc = cons_b.size();
    printf(" && nvars %d, n cv faces %d \n", n_vars, n_fc);
    assert(n_vars % 3 == 0);
    size_t n_verts = n_vars/3;

    std::vector<GRBVar> vars(n_vars);
    // initialize vars
    // GRBVar* vars = model.addVars(n_vars,, GRB_CONTINUOUS);
    for (size_t i = 0; i < n_vars; i++) {
        vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0., GRB_CONTINUOUS);
        vars[i].set(GRB_DoubleAttr_Start, x_0[i]); // this correct initialization?
    }

    // Integrate new variables
    model.update();

    //----------------------------------------------
    // 2. setup constraints
    //----------------------------------------------

    // coherent angles sum to π in each triangle
    for (size_t j = 0; j < n_fc; j++) {
        Eigen::VectorXd rowj = cons_A.row(j);
        double rhsj = cons_b[j];
        for (size_t i = 0; i < n_verts; i++){
            GRBLinExpr face_dist = rowj[0] * vars[3*i + 0] + 
                                   rowj[1] * vars[3*i + 1] + 
                                   rowj[2] * vars[3*i + 2];
            model.addConstr(face_dist <= rhsj);
        }
    }

    //----------------------------------------------
    // 3. setup energy
    //----------------------------------------------
    GRBQuadExpr objective;

    double constant_term = 0;
    // quad terms
    for (int k=0; k < Q.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(Q,k); it; ++it) {
            double val = it.value();
            size_t i = it.row(),
                   j = it.col();
            objective.addTerm(val, vars[i],  vars[j]);
        }
    }
    // linear terms
    for (int i = 0; i < n_vars; i++){
        objective.addTerm(g[i], vars[i]);
    }


    // set tolerance used for determining problem optimality
    // https://www.gurobi.com/documentation/current/refman/tolerances_and_user_scalin.html
    model.set(GRB_DoubleParam_OptimalityTol, 1e-8);

    // https://www.gurobi.com/documentation/current/refman/method.html
    // Algorithm used to solve continuous models or the initial root relaxation
    // of a MIP model. Options are:
    //   -1=automatic,
    //   0=primal simplex,
    //   1=dual simplex,
    //   2=barrier,
    //   3=concurrent,
    //   4=deterministic concurrent, and
    //   5=deterministic concurrent simplex.
    model.set(GRB_IntParam_Method, -1);
    model.setObjective(objective, GRB_MINIMIZE);
    // model.setObjectiveN(objective, 0, 1, 1); // main objective
    // model.set(GRB_IntParam_LogToConsole, 0);
    model.update();
    //----------------------------------------------
    // 4. solve problem
    //----------------------------------------------
    model.optimize();

    Eigen::VectorXd sol_x(x_0.size());
    for (int i = 0; i < n_vars; i++) {
        sol_x[i] = vars[i].get(GRB_DoubleAttr_X);
    }

    return sol_x;
}



///// MOSEK attempt
//     size_t n = x_0.size(),
//            m = cons_b.size();
//     Model::t M = new Model(); 
//     auto x = M->variable("x", n);
    
//     Matrix::t A_mat = Matrix::sparse(cons_A.)
//     // t >= |Ax-b|_2 where b is a parameter
    
//     auto b = M->parameter("b", m);

//   auto t = M->variable();
//   M->constraint(Expr::vstack(t, Expr::sub(Expr::mul(A, x), b)), Domain::inQCone());

//   // p_i >= |x_i|, i=1..n
//   auto p = M->variable(n);
//   M->constraint(Expr::hstack(p, x), Domain::inQCone());

//   // q >= |x|_2
//   auto q = M->variable();
//   M->constraint(Expr::vstack(q, x), Domain::inQCone());

//   // Objective, parametrized with lambda1, lambda2
//   // t + lambda1*sum(p) + lambda2*q
//   auto lambda1 = M->parameter("lambda1");
//   auto lambda2 = M->parameter("lambda2");
//   auto obj = Expr::add(new_array_ptr<Expression::t, 1>({t, Expr::mul(lambda1, Expr::sum(p)), Expr::mul(lambda2, q)}));
//   M->objective(ObjectiveSense::Minimize, obj);

//   // Return the ready model
//   return M;
// }
