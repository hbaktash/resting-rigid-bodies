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



#include <TinyAD/Support/GeometryCentral.hh>
#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>

#include <remesh_tools.h>

// #include <igl/arap.h>
#include <Eigen/Core>
#include "utils.h"
#include "geometry_utils.h"

#include "optimization.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;



// convertion stuff; Moved to utils

// tiny gradient stuff; Moved to utils


Vector3 barycentric(Vector3 p, Vector3 A, Vector3 B, Vector3 C);

class DeformationSolver{
    public:
       ManifoldSurfaceMesh *mesh;
       VertexPositionGeometry *old_geometry, *deformed_geometry;
       
       ManifoldSurfaceMesh *convex_mesh;
       VertexPositionGeometry *convex_geometry;

       // CP energy stuff
       size_t threshold_exceeded = 0;
       VertexData<SurfacePoint> closest_point_assignment; // hullMesh -> innerMesh
       VertexData<bool> vertex_only_assignment; // hullMesh -> bool
       VertexData<bool> marked_to_split; // hullMesh -> bool
       bool enforce_snapping = false;
       VertexData<Vertex> frozen_assignment; // hullMesh -> innerMesh
       Eigen::VectorX<bool> get_frozen_flags();
       Eigen::VectorXd get_frozen_x();

       VertexData<double> closest_point_distance; // hullMesh -> double
       VertexData<bool> CP_involvement;
       SparseMatrix<double> closest_point_operator;
       SparseMatrix<double> closest_point_flat_operator;
       size_t face_assignments = 0, edge_assignments = 0, vertex_assignments = 0;

       bool dynamic_remesh = true,
            curvature_weighted_CP = true,
            snap_after_split = true;

       double refinement_CP_threshold = 0.01,
              active_set_threshold = 0.08,
              split_robustness_threshold = 0.05;
       double init_bending_lambda = 1e0,
              final_bending_lambda = 1e3,
              init_membrane_lambda = 1e3,
              final_membrane_lambda = 1e0,
              init_CP_lambda = 1e1,
              final_CP_lambda = 1e9,
              reg_lambda = 0.;
       double final_barrier_lambda = 0.,
              init_barrier_lambda = 0.;
       // TODO: should this be different for every energy? 
       double internal_growth_p = 0.9, // p; for b-(a-b)p^{t}
              internal_pt = 1.;
       double get_scheduled_weight(double init_w, double final_w);

       int filling_max_iter = 50;
       // linear constraints
       DenseMatrix<double> constraint_matrix;
       Vector<double> constraint_rhs;

       // constructors
       DeformationSolver(ManifoldSurfaceMesh *old_mesh, VertexPositionGeometry *old_geometry,
                     ManifoldSurfaceMesh *convex_mesh, VertexPositionGeometry *convex_geometry);

       // bending energy
       double bending_energy(VertexPositionGeometry *new_geometry);
       // gradient of bending energy
       VertexData<Vector3> bending_energy_gradient(VertexPositionGeometry *new_geometry);
       // hessian of bending energy
       // DenseMatrix<double> bending_energy_hessian(VertexPositionGeometry *new_geometry);
       

       // closest point energy
       void assign_closest_vertices(VertexPositionGeometry *new_geometry, bool allow_multi_assignment = true);
       void assign_closest_points_barycentric(VertexPositionGeometry *new_geometry);
       // edge/face splits
       SurfacePoint get_robust_barycentric_point(SurfacePoint p, double threshold);
       bool split_barycentric_closest_points(VertexPositionGeometry *new_geometry);
       double closest_point_energy(VertexPositionGeometry *new_geometry);
       double closest_point_energy(Vector<double> flat_new_pos_mat);
       // gradient of CP energy
       DenseMatrix<double> closest_point_energy_gradient(VertexPositionGeometry *new_geometry);
       

       // barrier and CP computes
       Vector<double> get_CP_barrier_multiplier_vetor(double CP_involed_constant = 0.1);
       void build_constraint_matrix_and_rhs();
       std::tuple<double, DenseMatrix<double>, std::vector<DenseMatrix<double>>> get_log_barrier_stuff(DenseMatrix<double> new_pos_mat, Vector<double> weights);
       bool check_feasibility(DenseMatrix<double> new_pos_mat);
       double get_log_barrier_energy(DenseMatrix<double> new_pos_mat);
       // active set helpers
       DenseMatrix<bool> get_active_set_matrix(DenseMatrix<double> new_pos_mat, double active_threshold);
       double get_point_distance_to_convex_hull(Eigen::VectorXd p, bool from_faces = true);

       // rest geometry constants for elastic energies
       EdgeData<double> rest_dihedral_angles,
                        rest_bending_constant;
       FaceData<Eigen::Matrix2d> rest_membrane_I_inverted;
       FaceData<double> rest_face_areas;

       void update_bending_rest_constants();
       void update_membrane_rest_constants();
       
       auto get_tinyAD_bending_function(); // EdgeData<double> &rest_constants, EdgeData<double> &rest_dihedral_angles
       auto get_tinyAD_membrane_function(); // FaceData<Eigen::Matrix2d> &rest_tensors_inverted, FaceData<double> &rest_face_areas
       auto get_tinyAD_barrier_function();

       EdgeData<double> bending_per_edge;
       FaceData<double> membrane_per_vertex;
       // solver
       DenseMatrix<double> solve_for_bending(int visual_per_step = 0, 
                                             bool energy_plot = false, int* current_iter = nullptr, float** ys = nullptr);
       // void solve_qp(std::vector<Energy*> energies, std::vector<Constraint*> constraints);
       void print_energies_after_transform(Eigen::Matrix3d A);
       void test_my_barrier_vs_tinyAD();

       DenseMatrix<double> solve_for_center_of_mass(int visual_per_step = 0, 
                                                    bool energy_plot = false, int* current_iter = nullptr, float** ys = nullptr);
       Vector3 current_G, goal_G;
       double current_volume;

       VertexData<DenseMatrix<double>> per_vertex_G_jacobian(VertexPositionGeometry *tmp_geometry);
       Eigen::VectorXd flat_distance_multiplier(VertexPositionGeometry *tmp_geometry, bool from_faces);
       Eigen::MatrixXd per_vertex_G_derivative(VertexPositionGeometry *tmp_geometry);
       double init_G_lambda = 1e1,
              final_G_lambda = 1e3;
       DenseMatrix<double> solve_for_G(int visual_per_step = 0, 
                                       bool energy_plot = false, int* current_iter = nullptr, float** ys = nullptr);

};


void visualize_cv_cp_assignments();


// ARAP
// geometrycentral::DenseMatrix<double> get_ARAP_positions(DenseMatrix<double> old_pos_mat,
//                                        DenseMatrix<double> new_pos_mat,
//                                        DenseMatrix<double> init_sol,
//                                        ManifoldSurfaceMesh &inner_mesh,
//                                        Vector<int> hull_indices);

