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

#pragma once

// #include "monty.h"
// #include "fusion.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include <Eigen/Core>
// #include "gurobi_c++.h"
#include <osqp++.h>
// #include "monty.h"
// #include "fusion.h"

// using namespace mosek::fusion;
// using namespace monty;


#define ANSI_FG_MAGENTA "\x1b[35m"
#define ANSI_FG_YELLOW "\x1b[33m"
#define ANSI_FG_GREEN "\x1b[32m"
#define ANSI_FG_WHITE "\x1b[37m"
#define ANSI_FG_RED "\x1b[31m"
#define ANSI_RESET "\x1b[0m"


double line_search(Eigen::VectorXd x0, Eigen::VectorXd d, double f0, Eigen::VectorXd g, 
                            std::function<double(Eigen::VectorXd)> eval, double s_max = 1.0, 
                            double shrink = 0.8, int max_iters = 64);

Eigen::VectorXd solve_QP_with_ineq_GRB(Eigen::SparseMatrix<double> Q, Eigen::VectorXd g, Eigen::VectorXd x_0, 
                                       Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b, 
                                       Eigen::VectorX<bool> frozen_flags, Eigen::VectorXd frozen_x,
                                       Eigen::MatrixX<bool> active_set = Eigen::MatrixX<bool>::Zero(0,0));


Eigen::VectorXd solve_QP_with_ineq_OSQP(Eigen::SparseMatrix<double> Q, Eigen::VectorXd g, Eigen::VectorXd x_0, 
                                        Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b, 
                                        Eigen::VectorX<bool> frozen_flags, Eigen::VectorXd frozen_x,
                                        Eigen::MatrixX<bool> active_set = Eigen::MatrixX<bool>::Zero(0,0));

std::tuple<Eigen::SparseMatrix<double>, Eigen::VectorXd, Eigen::VectorXd> 
inequality_constraints_to_matrix(Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b, 
                                Eigen::VectorX<bool> frozen_flags, Eigen::VectorXd frozen_x,
                                Eigen::MatrixX<bool> active_set = Eigen::MatrixX<bool>::Zero(0,0));

// void build_QP_model_with_constraints(GRBModel &model, 
//                                      Eigen::VectorXd x_0, 
//                                      Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b, 
//                                      Eigen::VectorX<bool> frozen_flags, Eigen::VectorXd frozen_x,
//                                      Eigen::MatrixX<bool> active_set = Eigen::MatrixX<bool>::Zero(0,0));

// Eigen::VectorXd update_QP_objective_and_solve(GRBModel &model, 
//                                               Eigen::SparseMatrix<double> Q, Eigen::VectorXd g, Eigen::VectorXd x_0);


// Eigen::VectorXd solve_QP_with_ineq_MOSEK(Eigen::SparseMatrix<double> Q, Eigen::VectorXd g, Eigen::VectorXd x_0, 
//                                          Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b);
