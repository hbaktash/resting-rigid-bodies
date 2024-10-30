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

#include <Eigen/Core>

// using namespace geometrycentral;
// using namespace geometrycentral::surface;


Eigen::Vector3d to_eigen(const geometrycentral::Vector3& _v);
geometrycentral::Vector3 to_geometrycentral(const Eigen::Vector3d& _v);
geometrycentral::Vector<double> tinyAD_flatten(geometrycentral::DenseMatrix<double> mat);
geometrycentral::DenseMatrix<double> unflat_tinyAD(geometrycentral::Vector<double> flat_mat);
geometrycentral::SparseMatrix<double> tinyADify_barrier_hess(std::vector<geometrycentral::DenseMatrix<double>> hessians);
geometrycentral::Vector<double> vec32vec(geometrycentral::Vector3 v);
geometrycentral::Vector3 vec_to_GC_vec3(geometrycentral::Vector<double> vec);
geometrycentral::DenseMatrix<double> vertex_data_to_matrix(geometrycentral::surface::VertexData<geometrycentral::Vector3> positions);
geometrycentral::surface::VertexData<geometrycentral::Vector3> vertex_matrix_to_data(Eigen::MatrixXd positions, geometrycentral::surface::ManifoldSurfaceMesh& mesh);
geometrycentral::DenseMatrix<double> face_data_to_matrix(geometrycentral::surface::FaceData<geometrycentral::Vector3> fdata);
// Eigen::SparseMatrix<double> tinyADify_constraint_mat(Eigen::MatrixXd A);
double binomial_dist(int n, int k);
