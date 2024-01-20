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



#include <igl/arap.h>
#include <Eigen/Core>

using namespace geometrycentral;
using namespace geometrycentral::surface;



// convertion stuff

Vector<double> vec32vec(Vector3 v);
Vector3 vec_to_GC_vec3(Vector<double> vec);
DenseMatrix<double> vertex_data_to_matrix(VertexData<Vector3> positions);
DenseMatrix<double> face_data_to_matrix(FaceData<Vector3> fdata);

Vector<double> tinyAD_flatten(DenseMatrix<double> mat);
DenseMatrix<double> unflat_tinyAD(Vector<double> flat_mat);

SparseMatrix<double> tinyADify_barrier_hess(std::vector<DenseMatrix<double>> hessians);

// Gradient stuff
// from: https://www.sciencedirect.com/science/article/pii/S0167839607000891 
// Diherdral angle derivative; the angle <BAC on unit sphere
Vector3 dihedral_angle_grad_G(Vector3 G, Vector3 A, Vector3 B, Vector3 C);
Vector3 dihedral_angle_grad_A(Vector3 G, Vector3 A, Vector3 B, Vector3 C);
Vector3 dihedral_angle_grad_B(Vector3 G, Vector3 A, Vector3 B, Vector3 C);
Vector3 dihedral_angle_grad_C(Vector3 G, Vector3 A, Vector3 B, Vector3 C);

Vector3 barycentric(Vector3 p, Vector3 A, Vector3 B, Vector3 C);

class DeformationSolver{
    public:
        ManifoldSurfaceMesh *mesh;
        VertexPositionGeometry *old_geometry;
        
        ManifoldSurfaceMesh *convex_mesh;
        VertexPositionGeometry *convex_geometry;

        // CP energy stuff
        VertexData<SurfacePoint> closest_point_assignment;
        SparseMatrix<double> closest_point_operator;
        SparseMatrix<double> closest_point_flat_operator;

        bool one_time_CP_assignment = true;
        double CP_lambda = 10.0,
               CP_mu = 1.; // grow the CP lambda; since we want it to be zero in the end
        double barrier_init_lambda = 1e2,
               barrier_decay = 0.8;
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
        void assign_closest_points(VertexPositionGeometry *new_geometry);
        double closest_point_energy(VertexPositionGeometry *new_geometry);
        double closest_point_energy(Vector<double> flat_new_pos_mat);
        // gradient of CP energy
        DenseMatrix<double> closest_point_energy_gradient(VertexPositionGeometry *new_geometry);
        

        // constraints
        void build_constraint_matrix_and_rhs();
        std::tuple<double, DenseMatrix<double>, std::vector<DenseMatrix<double>>> get_log_barrier_stuff(DenseMatrix<double> new_pos_mat);

        // solver
        DenseMatrix<double> solve_for_bending(int visual_per_step = 0);
        // void solve_qp(std::vector<Energy*> energies, std::vector<Constraint*> constraints);
};



// ARAP
geometrycentral::DenseMatrix<double> get_ARAP_positions(DenseMatrix<double> old_pos_mat,
                                       DenseMatrix<double> new_pos_mat,
                                       DenseMatrix<double> init_sol,
                                       ManifoldSurfaceMesh &inner_mesh,
                                       Vector<int> hull_indices);

