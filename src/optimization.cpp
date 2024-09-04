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



Eigen::VectorXd line_search(Eigen::VectorXd x0, Eigen::VectorXd d, double f0, Eigen::VectorXd g, 
                            std::function<double(Eigen::VectorXd)> eval, double s_max, 
                            double shrink, int max_iters, double armijo_const){
    // Check input
    assert(x0.size() == g.size());
    if (s_max <= 0.0)
        throw std::runtime_error("Max step size not positive.");
    int i = 0;
    for (i = 0; i < max_iters; i++) {
        Eigen::VectorXd x_new = x0 + s_max * d;
        double f_new = eval(x_new);
        if (f_new <= f0 + armijo_const * s_max * d.dot(g)){
            std::cout << ANSI_FG_GREEN << " line search ended at iter " << i << "/" << max_iters << ANSI_RESET << std::endl;
            return x_new;
        }
        if (s_max > 1.0 && s_max * shrink < 1.0)
            s_max = 1.0;
        else
            s_max *= shrink;
    }
    std::cout << ANSI_FG_YELLOW << "Line search couldn't find improvement. Gradient max norm is" << g.cwiseAbs().maxCoeff() << ANSI_RESET << std::endl;
    return x0;
}


Eigen::VectorXd solve_QP_with_ineq_GRB(Eigen::SparseMatrix<double> Q, Eigen::VectorXd g, Eigen::VectorXd x_0, 
                                   Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b, 
                                   Eigen::VectorX<bool> frozen_flags, Eigen::VectorXd frozen_x,
                                   Eigen::MatrixX<bool> active_set){
                                        // Gurobi calling code borrowed from CoMISo
    // ******* Gurobi *******
    // build model
    // GRBEnv env     = GRBEnv();
    // GRBModel model = GRBModel(env);
    // build_QP_model_with_constraints(model, x_0, cons_A, cons_b, frozen_flags, frozen_x, active_set);
    // // set objective and solve
    // return update_QP_objective_and_solve(model, Q, g, x_0);
    return x_0;
}


Eigen::VectorXd solve_QP_with_ineq_OSQP(Eigen::SparseMatrix<double> Q, Eigen::VectorXd g, Eigen::VectorXd x_0, 
                                        Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b, 
                                        Eigen::VectorX<bool> frozen_indices, Eigen::VectorXd frozen_x,
                                        Eigen::MatrixX<bool> active_set){
    return x_0;
    // Eigen::SparseMatrix<double> objective_matrix(N, N);
    // Eigen::VectorXd objective_vector(N);
    // SparseMatrix<double> constraint_matrix(N + 2 * s2_mesh.nEdges(), N);

    // osqp::OsqpInstance instance;
    // instance.objective_matrix = objective_matrix;
    // instance.objective_vector.resize(2);
    // instance.objective_vector << 1.0, 0.0;
    // instance.constraint_matrix = constraint_matrix;
    // instance.lower_bounds.resize(1);
    // instance.lower_bounds << 1.0;
    // instance.upper_bounds.resize(1);
    // instance.upper_bounds << kInfinity;

    // osqp::OsqpSolver solver;
    // osqp::OsqpSettings settings;
    // // Edit settings if appropriate.
    // auto status = solver.Init(instance, settings);
    // // Assuming status.ok().
    // osqp::OsqpExitCode exit_code = solver.Solve();
    // // Assuming exit_code == OsqpExitCode::kOptimal.
    // double optimal_objective = solver.objective_value();
    // Eigen::VectorXd optimal_solution = solver.primal_solution();
}


// ======================= Gurobi =======================




// void build_QP_model_with_constraints(GRBModel &model, Eigen::VectorXd x_0, 
//                                      Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b, 
//                                      Eigen::VectorX<bool> frozen_flags, Eigen::VectorXd frozen_x,
//                                      Eigen::MatrixX<bool> active_set){
//     // log to console
//     model.set(GRB_IntParam_LogToConsole, 0);

//     //----------------------------------------------
//     // 0. set up gurobi environment
//     //----------------------------------------------

//     double time_limit = 1e10;
//     model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);

//     //----------------------------------------------
//     // 1. allocate variables
//     //----------------------------------------------
//     // GUROBI variables

//     size_t n_vars = x_0.size(),
//            n_fc = cons_b.size();
//     // printf(" && nvars %d, n cv faces %d \n", n_vars, n_fc);
//     assert(n_vars % 3 == 0);
//     size_t n_verts = n_vars/3;

//     std::vector<GRBVar> vars(n_vars);
//     // initialize vars
//     // GRBVar* vars = model.addVars(n_vars,, GRB_CONTINUOUS);
//     for (size_t i = 0; i < n_vars; i++) {
//         vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0., GRB_CONTINUOUS, "x" + std::to_string(i));
//         vars[i].set(GRB_DoubleAttr_Start, x_0[i]); // this correct initialization?
//     }

//     // Integrate new variables
//     model.update();

//     //----------------------------------------------
//     // 2. setup constraints
//     //----------------------------------------------

//     // add linear inequalities
//     for (size_t j = 0; j < n_fc; j++) {
//         Eigen::VectorXd rowj = cons_A.row(j);
//         double rhsj = cons_b[j];
//         for (size_t i = 0; i < n_verts; i++){
//             if (frozen_flags[i]) // skip frozen vertices
//                 continue; // add equality consts later
//             // GRBLinExpr face_dist = ;
//             if (active_set.size() > 1) // not default argument
//                 if (!active_set(j, i)) // skip inactive constraints
//                     continue;
//             model.addConstr(rowj[0] * vars[3*i + 0] + 
//                             rowj[1] * vars[3*i + 1] + 
//                             rowj[2] * vars[3*i + 2] <= rhsj);
//         }
//     }
//     // add frozen inequalities
//     for (size_t i = 0; i < n_verts; i++){
//         if (frozen_flags[i])
//             model.addConstr(vars[i] == frozen_x[i]);
//     }

//     model.update();
// }


// Eigen::VectorXd update_QP_objective_and_solve(GRBModel &model, 
//                                               Eigen::SparseMatrix<double> Q, Eigen::VectorXd g, Eigen::VectorXd x_0){
//     //----------------------------------------------
//     // retreive variables
//     //----------------------------------------------
//     size_t n_vars = x_0.size();
//     std::vector<GRBVar> vars(n_vars);
//     for (int i = 0; i < n_vars; i++) {
//         vars[i] = model.getVarByName("x" + std::to_string(i));
//         vars[i].set(GRB_DoubleAttr_Start, x_0[i]);
//     }

//     //----------------------------------------------
//     // 3. setup energy
//     //----------------------------------------------
//     GRBQuadExpr objective;

//     double constant_term = 0;
//     // quad terms
//     for (int k=0; k < Q.outerSize(); ++k){
//         for (Eigen::SparseMatrix<double>::InnerIterator it(Q,k); it; ++it) {
//             double val = it.value();
//             size_t i = it.row(),
//                    j = it.col();
//             objective.addTerm(val, vars[i], vars[j]);
//         }
//     }
//     // linear terms
//     for (int i = 0; i < n_vars; i++){
//         objective.addTerm(g[i], vars[i]);
//     }


//     // set tolerance used for determining problem optimality
//     // https://www.gurobi.com/documentation/current/refman/tolerances_and_user_scalin.html
//     model.set(GRB_DoubleParam_OptimalityTol, 1e-8);

//     // https://www.gurobi.com/documentation/current/refman/method.html
//     // Algorithm used to solve continuous models or the initial root relaxation
//     // of a MIP model. Options are:
//     //   -1=automatic,
//     //   0=primal simplex,
//     //   1=dual simplex,
//     //   2=barrier,
//     //   3=concurrent,
//     //   4=deterministic concurrent, and
//     //   5=deterministic concurrent simplex.
//     model.set(GRB_IntParam_Method, -1);
//     // model.set(GRB_IntParam_BarOrder, 1);
//     model.setObjective(objective, GRB_MINIMIZE);
//     // model.setObjectiveN(objective, 0, 1, 1); // main objective
//     model.update();
//     //----------------------------------------------
//     // 4. solve problem
//     //----------------------------------------------
//     try {
//         model.optimize();
//         int optimstatus = model.get(GRB_IntAttr_Status);
//         if (optimstatus == GRB_OPTIMAL) { // Solved!
//             std::cout << ANSI_FG_GREEN << "Optimal objective: " << model.get(GRB_DoubleAttr_ObjVal) << ANSI_RESET << std::endl;
//             Eigen::VectorXd sol_x(x_0.size());
//             for (int i = 0; i < n_vars; i++) {
//                 sol_x[i] = vars[i].get(GRB_DoubleAttr_X);
//             }
//             return sol_x;
//         // Not solved
//         } else if (optimstatus == GRB_INFEASIBLE) { 
//             std::cout << ANSI_FG_RED << "Model is infeasible" << ANSI_RESET << std::endl;
//         } else if (optimstatus == GRB_UNBOUNDED) {
//             std::cout << ANSI_FG_RED << "Model is unbounded" << ANSI_RESET << std::endl;
//         } else {
//             std::cout << ANSI_FG_RED << "Optimization was stopped with status = "<< ANSI_RESET << optimstatus << std::endl;
//         }
//     } catch (GRBException e) {
//         std::cout << ANSI_FG_RED << "Error code = " << e.getErrorCode() << ANSI_RESET << std::endl;
//         std::cout << ANSI_FG_YELLOW << e.getMessage() << ANSI_RESET << std::endl;
//     } catch (...) {
//         std::cout << ANSI_FG_RED << "Error during optimization" << ANSI_RESET << std::endl;
//     }
//     return x_0;                         
// }




// ======================= MOSEK =======================


// Eigen::VectorXd solve_QP_with_ineq_MOSEK(Eigen::SparseMatrix<double> Q, Eigen::VectorXd g, Eigen::VectorXd x_0, 
//                                          Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b){
// //          min_x  |Ax-b|_2 + lambda1*|x|_1 + lambda2*|x|_2

//     size_t n = x_0.size(),
//             m = cons_b.size();
//     Model::t M = new Model(); 

//     // add variables
//     Variable::t x = M->variable("x", n);

//     // add constraints
//     auto b = M->parameter("b", m); // cons_b
//     M->constraint(Expr::vstack(t, Expr::sub(Expr::mul(cons_A, x), b)), Domain::inQCone());

//     // p_i >= |x_i|, i=1..n
//     auto p = M->variable(n);
//     M->constraint(Expr::hstack(p, x), Domain::inQCone());

//     // q >= |x|_2
//     auto q = M->variable();
//     M->constraint(Expr::vstack(q, x), Domain::inQCone());

//     // Objective, parametrized with lambda1, lambda2
//     // t + lambda1*sum(p) + lambda2*q
//     auto lambda1 = M->parameter("lambda1");
//     auto lambda2 = M->parameter("lambda2");
//     auto obj = Expr::add(new_array_ptr<Expression::t, 1>({t, Expr::mul(lambda1, Expr::sum(p)), Expr::mul(lambda2, q)}));
//     M->objective(ObjectiveSense::Minimize, obj);

//     // Return the ready model
// }

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
