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

#include <igl/arap.h>
#include <Eigen/Core>

using namespace geometrycentral;
using namespace geometrycentral::surface;


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

        VertexData<SurfacePoint> closest_point_assignment;
        SparseMatrix<double> closest_point_operator;
        // constructors
        DeformationSolver(ManifoldSurfaceMesh *old_mesh, VertexPositionGeometry *old_geometry,
                          ManifoldSurfaceMesh *convex_mesh, VertexPositionGeometry *convex_geometry);

        // bending energy
        double bending_energy(VertexPositionGeometry *new_geometry);
        // gradient of energy X
        VertexData<Vector3> bending_energy_gradient(VertexPositionGeometry *new_geometry);
        // hessian of energy X
        DenseMatrix<double> bending_energy_hessian(VertexPositionGeometry *new_geometry);
        
        // closest point energy
        void assign_closest_points(VertexPositionGeometry *new_geometry);
        double closest_point_energy(VertexPositionGeometry *new_geometry);
        // gradient of CP energy
        VertexData<Vector3> closest_point_energy_gradient(VertexPositionGeometry *new_geometry);
        

        // constraints

        // QP solver for a set of energies
        // void solve_qp(std::vector<Energy*> energies, std::vector<Constraint*> constraints);
};



// ARAP
geometrycentral::DenseMatrix<double> get_ARAP_positions(DenseMatrix<double> old_pos_mat,
                                       DenseMatrix<double> new_pos_mat,
                                       DenseMatrix<double> init_sol,
                                       ManifoldSurfaceMesh &inner_mesh,
                                       Vector<int> hull_indices);

