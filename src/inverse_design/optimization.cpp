#include "optimization.h"



double line_search(Eigen::VectorXd x0, Eigen::VectorXd d, double f0, Eigen::VectorXd g, 
                            std::function<double(Eigen::VectorXd)> eval, double s_max, 
                            double shrink, int max_iters){
    // Check input
    assert(x0.size() == g.size());
    if (s_max < 0.0)
        throw std::runtime_error("Max step size not positive.");
    int i = 0;
    for (i = 0; i < max_iters; i++) {
        Eigen::VectorXd x_new = x0 + s_max * d;
        double f_new = eval(x_new);
        if (f_new <= f0){
            std::cout << ANSI_FG_GREEN << " line search ended at iter " << i << "/" << max_iters << " with s = " << s_max << ANSI_RESET << std::endl;
            return s_max;
        }
        if (s_max > 1.0 && s_max * shrink < 1.0)
            s_max = 1.0;
        else
            s_max *= shrink;
    }
    std::cout << ANSI_FG_YELLOW << "Line search couldn't find improvement. Gradient max norm is" << g.cwiseAbs().maxCoeff() << ANSI_RESET << std::endl;
    return 0.;
}


Eigen::VectorXd solve_QP_with_ineq_GRB(Eigen::SparseMatrix<double> Q, Eigen::VectorXd g, Eigen::VectorXd x_0, 
                                   Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b, 
                                   Eigen::VectorX<bool> frozen_flags, Eigen::VectorXd frozen_x,
                                   Eigen::MatrixX<bool> active_set){
    // ******* Gurobi *******
    // build model
    GRBEnv env     = GRBEnv();
    GRBModel model = GRBModel(env);
    build_GRB_QP_model_with_constraints(model, x_0, cons_A, cons_b, frozen_flags, frozen_x, active_set);
    // set objective and solve
    return update_GRB_QP_objective_and_solve(model, Q, g, x_0);
    // return x_0;
}

std::tuple<Eigen::SparseMatrix<double>, Eigen::VectorXd, Eigen::VectorXd> 
inequality_constraints_to_matrix(Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b, 
                                Eigen::VectorX<bool> frozen_flags, Eigen::VectorXd frozen_x,
                                Eigen::MatrixX<bool> active_set){
    const double kInfinity = std::numeric_limits<double>::infinity();
    
    // build matrix
    size_t n_vars = frozen_x.size(),
           n_fc = cons_b.size();
    size_t n_verts = n_vars/3;
    // add linear inequalities
    std::vector<double> lb_list, ub_list;
    std::vector<Eigen::Triplet<double>> tripletList;

    size_t const_counter = 0;
    for (size_t j = 0; j < n_fc; j++) {
        Eigen::VectorXd rowj = cons_A.row(j);
        double rhsj = cons_b[j];
        for (size_t i = 0; i < n_verts; i++){
            if (frozen_flags[i]) // skip frozen vertices
                continue; // add equality consts later
            // GRBLinExpr face_dist = ;
            if (active_set.size() > 1) // not default argument
                if (!active_set(j, i)) // skip inactive constraints
                    continue;
            tripletList.push_back(Eigen::Triplet<double>(const_counter, 3*i + 0, rowj[0]));
            tripletList.push_back(Eigen::Triplet<double>(const_counter, 3*i + 1, rowj[1]));
            tripletList.push_back(Eigen::Triplet<double>(const_counter, 3*i + 2, rowj[2]));
            lb_list.push_back(-kInfinity);
            ub_list.push_back(rhsj);
            const_counter++;
            // model.addConstr(rowj[0] * vars[3*i + 0] + 
            //                 rowj[1] * vars[3*i + 1] + 
            //                 rowj[2] * vars[3*i + 2] <= rhsj);
        }
    }
    // add frozen inequalities
    for (size_t i = 0; i < n_vars; i++){
        if (frozen_flags[i]){
            tripletList.push_back(Eigen::Triplet<double>(const_counter, i, 1.)); //tODO
            lb_list.push_back(frozen_x[i]);
            ub_list.push_back(frozen_x[i]);
            const_counter++;
            // model.addConstr(vars[i] == frozen_x[i]);
        }
    }
    Eigen::SparseMatrix<double> constraint_matrix(const_counter, n_vars);
    constraint_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
    Eigen::VectorXd lower_bounds = Eigen::Map<Eigen::VectorXd>(lb_list.data(), lb_list.size());
    Eigen::VectorXd upper_bounds = Eigen::Map<Eigen::VectorXd>(ub_list.data(), ub_list.size());
    return std::make_tuple(constraint_matrix, lower_bounds, upper_bounds);
}

Eigen::VectorXd solve_QP_with_ineq_OSQP(Eigen::SparseMatrix<double> Q, Eigen::VectorXd g, Eigen::VectorXd x_0, 
                                        Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b, 
                                        Eigen::VectorX<bool> frozen_flags, Eigen::VectorXd frozen_x,
                                        Eigen::MatrixX<bool> active_set){
    // generic LP solve
    // Eigen::VectorXd build_and_solve_LP(Eigen::SparseMatrix<double> objective_matrix, Eigen::VectorXd 
    //                                 objective_vector, Eigen::SparseMatrix<double> constraint_matrix, 
    //                                 Eigen::VectorXd lower_bounds, Eigen::VectorXd upper_bounds,
    //                                 int max_iter, double err_tol){
        // build OSQP model
    assert(lower_bounds.size() == constraint_matrix.rows());
    assert(upper_bounds.size() == constraint_matrix.rows());

    osqp::OsqpInstance instance;
    instance.objective_matrix = Q;
    instance.objective_vector = g;

    Eigen::SparseMatrix<double> constraint_matrix;
    Eigen::VectorXd lower_bounds, upper_bounds;
    std::tie(constraint_matrix, lower_bounds, upper_bounds) = inequality_constraints_to_matrix(cons_A, cons_b, frozen_flags, frozen_x, active_set);
    instance.constraint_matrix = constraint_matrix;
    instance.lower_bounds = lower_bounds;
    instance.upper_bounds = upper_bounds;
    
    osqp::OsqpSolver osqp_solver;
    osqp::OsqpSettings osqp_settings; 
    osqp_settings.verbose = false;
    

    // Edit settings
    osqp_settings.max_iter = 5e4;
    auto status = osqp_solver.Init(instance, osqp_settings);
    // Assuming status.ok().
    osqp::OsqpExitCode exit_code = osqp_solver.Solve();
    // Assuming exit_code == OsqpExitCode::kOptimal.
    std::cout << ANSI_FG_MAGENTA << " \t QP solve done,  Exit code: " << ToString(exit_code) << 
                                    " \n\t total steps: " << osqp_solver.iterations() << ANSI_RESET << std::endl;
    double optimal_objective = osqp_solver.objective_value();
    Eigen::VectorXd optimal_solution = osqp_solver.primal_solution();
    return optimal_solution;
}






// ======================= Gurobi =======================




void build_GRB_QP_model_with_constraints(GRBModel &model, Eigen::VectorXd x_0, 
                                     Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b, 
                                     Eigen::VectorX<bool> frozen_flags, Eigen::VectorXd frozen_x,
                                     Eigen::MatrixX<bool> active_set){
    // log to console
    model.set(GRB_IntParam_LogToConsole, 0);

    //----------------------------------------------
    // 0. set up gurobi environment
    //----------------------------------------------

    double time_limit = 1e10;
    model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);

    //----------------------------------------------
    // 1. allocate variables
    //----------------------------------------------
    // GUROBI variables

    size_t n_vars = x_0.size(),
           n_fc = cons_b.size();
    // printf(" && nvars %d, n cv faces %d \n", n_vars, n_fc);
    assert(n_vars % 3 == 0);
    size_t n_verts = n_vars/3;

    std::vector<GRBVar> vars(n_vars);
    // initialize vars
    // GRBVar* vars = model.addVars(n_vars,, GRB_CONTINUOUS);
    for (size_t i = 0; i < n_vars; i++) {
        vars[i] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0., GRB_CONTINUOUS, "x" + std::to_string(i));
        vars[i].set(GRB_DoubleAttr_Start, x_0[i]); // this correct initialization?
    }

    // Integrate new variables
    model.update();

    //----------------------------------------------
    // 2. setup constraints
    //----------------------------------------------

    // add linear inequalities
    for (size_t j = 0; j < n_fc; j++) {
        Eigen::VectorXd rowj = cons_A.row(j);
        double rhsj = cons_b[j];
        for (size_t i = 0; i < n_verts; i++){
            if (frozen_flags[i]) // skip frozen vertices
                continue; // add equality consts later
            // GRBLinExpr face_dist = ;
            if (active_set.size() > 1) // not default argument
                if (!active_set(j, i)) // skip inactive constraints
                    continue;
            model.addConstr(rowj[0] * vars[3*i + 0] + 
                            rowj[1] * vars[3*i + 1] + 
                            rowj[2] * vars[3*i + 2] <= rhsj);
        }
    }
    // add frozen inequalities
    for (size_t i = 0; i < n_vars; i++){
        if (frozen_flags[i])
            model.addConstr(vars[i] == frozen_x[i]);
    }

    model.update();
}


Eigen::VectorXd update_GRB_QP_objective_and_solve(GRBModel &model, 
                                              Eigen::SparseMatrix<double> Q, Eigen::VectorXd g, Eigen::VectorXd x_0){
    //----------------------------------------------
    // retreive variables
    //----------------------------------------------
    size_t n_vars = x_0.size();
    std::vector<GRBVar> vars(n_vars);
    for (int i = 0; i < n_vars; i++) {
        vars[i] = model.getVarByName("x" + std::to_string(i));
        vars[i].set(GRB_DoubleAttr_Start, x_0[i]);
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
            objective.addTerm(val, vars[i], vars[j]);
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
    // model.set(GRB_IntParam_BarOrder, 1);
    model.setObjective(objective, GRB_MINIMIZE);
    // model.setObjectiveN(objective, 0, 1, 1); // main objective
    model.update();
    //----------------------------------------------
    // 4. solve problem
    //----------------------------------------------
    try {
        model.optimize();
        int optimstatus = model.get(GRB_IntAttr_Status);
        if (optimstatus == GRB_OPTIMAL) { // Solved!
            std::cout << ANSI_FG_GREEN << "Optimal objective: " << model.get(GRB_DoubleAttr_ObjVal) << ANSI_RESET << std::endl;
            Eigen::VectorXd sol_x(x_0.size());
            for (int i = 0; i < n_vars; i++) {
                sol_x[i] = vars[i].get(GRB_DoubleAttr_X);
            }
            return sol_x;
        // Not solved
        } else if (optimstatus == GRB_INFEASIBLE) { 
            std::cout << ANSI_FG_RED << "Model is infeasible" << ANSI_RESET << std::endl;
        } else if (optimstatus == GRB_UNBOUNDED) {
            std::cout << ANSI_FG_RED << "Model is unbounded" << ANSI_RESET << std::endl;
        } else {
            std::cout << ANSI_FG_RED << "Optimization was stopped with status = "<< ANSI_RESET << optimstatus << std::endl;
        }
    } catch (GRBException e) {
        std::cout << ANSI_FG_RED << "Error code = " << e.getErrorCode() << ANSI_RESET << std::endl;
        std::cout << ANSI_FG_YELLOW << e.getMessage() << ANSI_RESET << std::endl;
    } catch (...) {
        std::cout << ANSI_FG_RED << "Error during optimization" << ANSI_RESET << std::endl;
    }
    return x_0;
}




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




// Eigen::VectorXd 
// solve_with_trust_region_Newton(Eigen::SparseMatrix<double> Q, Eigen::VectorXd g, Eigen::VectorXd x_0, 
//                                Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b, 
//                                std::vector<int> &mask, // mask for boundary constraints
//                                double initRadius,
//                                double maxRadius,
//                                double stopEpsilon = 1e-8,
//                                const int maxIterations = 100,
//                                const int cgIterations = 100,
//                                double eta = 0.25,
//                                bool quiet = false ){
//     // Trust region Newton method
//     // from Nocedal & Wright book: numerical optimization
//     // Algorithm 4.1
    
//     const std::vector<int> *_bdryMask;


//     double trRadius = initRadius;
//     double eta_k, eps_k, rho_k;

//     int n = x_0.size();

//     Eigen::VectorXd x_k = x_0;

//     Eigen::VectorXd p_k( x_0.size());
//     p_k.setZero();

//     Eigen::VectorXd tmp_x_k( x_0.size());
//     double F_k, tmp_F_k;
//     double m_red, f_red; // predicted and actual reduction

//     Eigen::VectorXd grad_F_k;
//     Eigen::MatrixXd Hess_F_k;
//     // _F.apply( x_k, F_k );
//     // _DF.apply( x_k, grad_F_k );
//     // _D2F.apply( x_k, Hess_F_k );
//     if ( _bdryMask )
//       applyMaskToSymmetricMatrixAndVector<Eigen::MatrixXd, Eigen::VectorXd>( *_bdryMask, Hess_F_k, grad_F_k );

//     double initGradNorm = grad_F_k.norm();

//     Eigen::VectorXd diagonal( x_0.size());
//     Eigen::MatrixXd Dinv( x_0.size(), x_0.size());
//     for ( int i = 0; i < n; i++ ) {
//       Dinv.coeffRef( i, i ) = 1;
//     }

//     Eigen::VectorXd c( n );
//     Eigen::MatrixXd H( n, n );

//     auto t_start_eval = std::chrono::high_resolution_clock::now();
//     auto t_end_eval = std::chrono::high_resolution_clock::now();


//     for ( int k = 0; k < maxIterations; k++ ) {
//       // Step 1: Solve trust-region subproblem
      
//       bool _diagonalPreconditioning = true; // TODO figure out
//       double eps_red = 100 * std::numeric_limits<double>::epsilon();
//       if ( _diagonalPreconditioning ) {
//         diagonal = Hess_F_k.diagonal();
//         for ( int i = 0; i < n; i++ ) {
//           if ( std::abs( diagonal[i] ) > eps_red )
//             Dinv.coeffRef( i, i ) = 1 / diagonal[i];
//           else
//             Dinv.coeffRef( i, i ) = 0;
//         }
//       }

//       c.noalias() = Dinv * grad_F_k;
//       H = Dinv * Hess_F_k * Dinv;

//       // Compute forcing sequence / epsilon
//       eta_k = std::min( 0.5, std::sqrt( c.norm()));
//       eps_k = eta_k * c.norm();

//       // Apply trust-region subproblem solver
//       SolverStatus<ConfiguratorType> trsolverStatus = solveTrustRegionSubproblem( H, c, trRadius, eps_k, p_k );
//       auto t_end_inner = std::chrono::high_resolution_clock::now();
//       this->status.additionalTimings["Subproblem"] += std::chrono::duration<double, std::milli>(
//               t_end_inner - t_start_inner ).count();

//       double pkn = p_k.norm();

//       p_k = Dinv * p_k;

//       // Step 2: Determine reduction ration
//       tmp_x_k = x_k + p_k; // temporary new iterate
//       _F.apply( tmp_x_k, tmp_F_k ); // temporary new function value

//       m_red = -grad_F_k.dot( p_k ) - 0.5 * p_k.dot( Hess_F_k * p_k ); // predicted reduction i.e. in the quadratic model
//       f_red = (F_k - tmp_F_k);

//       if ((std::abs( f_red ) < eps_red && std::abs( m_red ) < eps_red) || std::abs( f_red - m_red ) < eps_red ) {
//         if ( !_quiet )
//           std::cout << " -- TRN -- Iter " << std::setw( 3 ) << k << ": " << "Cutoff active in rho-computation"
//                     << std::endl;
//         rho_k = 1.;
//       }
//       else
//         rho_k = f_red / m_red; // actual over predicted reduction

//       // Step 3: Update trust region radius
//       if ( rho_k < 0.25 )
//         trRadius = trRadius / 4.;
//       else if ( rho_k > 0.75 && std::abs( pkn - trRadius ) <= eps_red )
//         trRadius = std::min( 2 * trRadius, _maxRadius );

//       // Step 4: Accept or decline new iterate
//       if ( rho_k > _eta ) {
//         x_k = tmp_x_k;
//         F_k = tmp_F_k;

//         t_start_eval = std::chrono::high_resolution_clock::now();
//         _DF.apply( x_k, grad_F_k );
//         _D2F.apply( x_k, Hess_F_k );
//         if ( _bdryMask )
//           applyMaskToSymmetricMatrixAndVector<Eigen::MatrixXd, Eigen::VectorXd>( *_bdryMask, Hess_F_k, grad_F_k );
//         t_end_eval = std::chrono::high_resolution_clock::now();
//         this->status.additionalTimings["Evaluation"] += std::chrono::duration<double, std::milli>(
//                 t_end_eval - t_start_eval ).count();

//         if ( grad_F_k.norm() < _stopEpsilon ) {
//           if ( !_quiet )
//             std::cout << " -- TRN -- Iter " << std::setw( 3 ) << k << ": " << "Gradient norm below epsilon."
//                     << std::endl;
//           auto t_end = std::chrono::high_resolution_clock::now();
//           this->status.totalTime += std::chrono::duration<double, std::milli>( t_end - t_start ).count();
//           break;
//         }
//       }

//       if ( trRadius < _minRadius ) {
//         if ( !_quiet )
//           std::cout << " -- TRN -- Iter " << std::setw( 3 ) << k << ": " << "Trust region too small." << std::endl;
//         auto t_end = std::chrono::high_resolution_clock::now();
//         this->status.totalTime += std::chrono::duration<double, std::milli>( t_end - t_start ).count();
//         break;
//       }

//       if ( p_k.template lpNorm<Eigen::Infinity>() < _minStepsize ) {
//         if ( !_quiet )
//           std::cout << " -- TRN -- Iter " << std::setw( 3 ) << k << ": " << "Step size too small." << std::endl;
//         auto t_end = std::chrono::high_resolution_clock::now();
//         this->status.totalTime += std::chrono::duration<double, std::milli>( t_end - t_start ).count();
//         break;
//       }

//     }
// }