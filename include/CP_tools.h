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
#include "geometrycentral/surface/surface_point.h"


// #include <igl/arap.h>
#include <Eigen/Core>
#include "utils.h"
#include "geometry_utils.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;


// static version of CP operators and energy and assignments
double get_raw_CP_energy(ManifoldSurfaceMesh &concave_mesh, Eigen::VectorXd concave_poses_tiny_AD_flattened, 
                       ManifoldSurfaceMesh &convex_mesh, VertexPositionGeometry &convex_geometry);
// std::pair<Eigen::MatrixXd, Eigen::VectorXd> 
// get_hull_constraint_matrix_and_rhs(ManifoldSurfaceMesh* convex_mesh, VertexPositionGeometry* convex_geometry){
//     size_t nf = convex_mesh->nFaces();
//     Eigen::MatrixXd constraint_matrix = Eigen::MatrixXd::Zero(nf, 3);
//     Eigen::VectorXd constraint_rhs = Eigen::VectorXd::Zero(nf);
//     for (Face f: convex_mesh->faces()){
//         Vector3 face_normal = convex_geometry->faceNormal(f); // assume outward normals
//         constraint_matrix.row(f.getIndex()) = vec32vec(face_normal);
//         Vector3 point_on_face_normal = convex_geometry->inputVertexPositions[f.halfedge().vertex()];
//         constraint_rhs(f.getIndex()) = dot(face_normal, point_on_face_normal);
//     }
//     return {constraint_matrix, constraint_rhs};
// }


// double 
// point_distance_to_convex_faces(Eigen::VectorXd p, 
//                                Eigen::MatrixXd constraint_matrix, Eigen::VectorXd constraint_rhs){
//     double distance = 0.;
//     Eigen::VectorXd Np = constraint_matrix * p;
//     Eigen::MatrixXd diff_ij = constraint_rhs - Np; // col j is N*P_j - rhs; i.e. -f_i(P_j) which should be positive
//     distance = diff_ij.minCoeff();
//     return distance;
// }


