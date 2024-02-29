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
#include "gurobi_c++.h"

// using namespace mosek::fusion;
// using namespace monty;

Eigen::VectorXd solve_QP_with_ineq(Eigen::SparseMatrix<double> Q, Eigen::VectorXd g, Eigen::VectorXd x_0, 
                                   Eigen::MatrixXd cons_A, Eigen::VectorXd cons_b);