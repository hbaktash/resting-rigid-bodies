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

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/numerical/linear_solvers.h"

#include <igl/arap.h>
#include <Eigen/Core>

using namespace geometrycentral;
using namespace geometrycentral::surface;


geometrycentral::DenseMatrix<double> get_ARAP_positions(geometrycentral::DenseMatrix<double> old_pos_mat,
                                       geometrycentral::DenseMatrix<double> new_pos_mat,
                                       geometrycentral::DenseMatrix<double> init_sol,
                                       ManifoldSurfaceMesh &inner_mesh,
                                       geometrycentral::Vector<int> hull_indices);

SparseMatrix<double> foo();

