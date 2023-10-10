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

// #include "geometrycentral/surface/manifold_surface_mesh.h"
// #include "geometrycentral/surface/vertex_position_geometry.h"
// #include "geometrycentral/surface/surface_point.h"
#include "boundary_tools.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;




// Gradient stuff
// from: https://www.sciencedirect.com/science/article/pii/S0167839607000891 
// Diherdral angle derivative; the angle <BAC on unit sphere
Vector3 dihedral_angle_grad_G(Vector3 G, Vector3 A, Vector3 B, Vector3 C);
Vector3 dihedral_angle_grad_A(Vector3 G, Vector3 A, Vector3 B, Vector3 C);
Vector3 dihedral_angle_grad_B(Vector3 G, Vector3 A, Vector3 B, Vector3 C);
Vector3 dihedral_angle_grad_C(Vector3 G, Vector3 A, Vector3 B, Vector3 C);


class InverseSolver{
    public:
        Forward3DSolver* forwardSolver;
        BoundaryBuilder* boundaryBuilder;

        FaceData<double> goal_area;

        InverseSolver(){}
        InverseSolver(BoundaryBuilder* boundaryBuilder);

        // distribution goals
        void set_fair_distribution();

        // gradient computation; assuming regularity
        // G grad; vertices frozen
        FaceData<Vector3> per_face_G_gradient;
        void find_per_face_G_grads(bool check_FD = false);
        Vector3 find_total_g_grad();
        // vertex grad; G frozen
        FaceData<VertexData<Vector3>> per_face_per_vertex_gradient;
        void find_per_face_per_vertex_grads(bool check_FD = false);
        VertexData<Vector3> find_per_vertex_total_grads();
        // G is dependent of Geometry; uniform mass
        FaceData<VertexData<Vector3>> per_face_per_vertex_uni_mass_gradient;
        void find_per_face_per_vertex_uni_mass_grads();
        VertexData<Vector3> find_per_vertex_uni_mass_total_grads();
        
};