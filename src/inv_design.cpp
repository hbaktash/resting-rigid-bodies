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

#include "inv_design.h"


// Gradient stuff
// formula source in header file
Vector3 dihedral_angle_grad_G(Vector3 G, Vector3 A, Vector3 B, Vector3 C){
    // A is in the middle (on the unit sphere)
    Vector3 GA = A - G,
            GB = B - G,
            GC = C - G;
    // un-normalized normals
    // order: B A C
    Vector3 nB_un = cross(GA, B - A),
            nC_un = cross(GA, A - C);
    Vector3 grad1 = nB_un * dot(-GA, B - A)/nB_un.norm2(),
            grad2 = nC_un * dot(GA, A - C)/nC_un.norm2();
    return -(grad1 + grad2)/GA.norm();
}
// same as G, replacing G and A
Vector3 dihedral_angle_grad_A(Vector3 G, Vector3 A, Vector3 B, Vector3 C){
    return dihedral_angle_grad_G(A, G, B, C);
}
// wing side vertex B; kinda simple
Vector3 dihedral_angle_grad_B(Vector3 G, Vector3 A, Vector3 B, Vector3 C){
    Vector3 nB_un = cross(A - G, B - A);
    return 0.5 * nB_un * (G-A).norm()/nB_un.norm2();
}
// wing side vertex C; kinda simple
Vector3 dihedral_angle_grad_C(Vector3 G, Vector3 A, Vector3 B, Vector3 C){
    Vector3 nC_un = cross(A - G, A - C);
    return 0.5 * nC_un * (G-A).norm()/nC_un.norm2();
}

// optimization stuff

InverseSolver::InverseSolver(BoundaryBuilder* boundaryBuilder){
    this->boundaryBuilder = boundaryBuilder;
    forwardSolver = boundaryBuilder->forward_solver;
    set_fair_distribution();
}

void InverseSolver::set_fair_distribution() {
    goal_area = FaceData<double>(*forwardSolver->hullMesh, 
                                         1./(double)forwardSolver->hullMesh->nFaces());
}

void InverseSolver::find_per_face_G_grads(bool check_FD){
    Vector3 zero_vec = Vector3::zero();
    Vector3 G = forwardSolver->get_G();
    per_face_G_gradient = FaceData<Vector3>(*forwardSolver->hullMesh, zero_vec);
    for (Face f: forwardSolver->hullMesh->faces()){
        Vector3 f_g_grad = zero_vec;
        // assuming "regularity"
        for (Halfedge he: f.adjacentHalfedges()){
            Vertex v0 = he.tailVertex(),
                   v1 = he.tipVertex(),
                   v2 = he.next().tipVertex();
            Vector3 p0 = forwardSolver->hullGeometry->inputVertexPositions[v0], 
                    p1 = forwardSolver->hullGeometry->inputVertexPositions[v1], 
                    p2 = forwardSolver->hullGeometry->inputVertexPositions[v2];
            Vector3 tmp_angle_G_grad = dihedral_angle_grad_G(G, p1, p0, p2);
            f_g_grad += tmp_angle_G_grad;
        }
        per_face_G_gradient[f] = f_g_grad;
    }
    if (check_FD){
        printf("FD check: \n");
        double step = 1e-6;
        Vector3 e_x({1., 0. ,0.}),
                e_y({0., 1. ,0.}),
                e_z({0., 0. ,1.});
        Vector<double> old_face_areas = boundaryBuilder->face_region_area.toVector();
        for (Vector3 dG: {e_x, e_y, e_z}){
            std::cout << "dG: "<< dG <<"\n";
            Vector3 tmp_G = G + dG * step;
            forwardSolver->set_G(tmp_G);
            forwardSolver->initialize_pre_computes();
            // boundaryBuilder->forward_solver = forwardSolver;
            boundaryBuilder->build_boundary_normals();
            Vector<double> new_face_areas = boundaryBuilder->face_region_area.toVector();
            for (Face f: forwardSolver->hullMesh->faces()){
                double proj_grad = dot(per_face_G_gradient[f], dG),
                       fd_grad = (new_face_areas[f.getIndex()] - old_face_areas[f.getIndex()])/step;
                printf(" proj grad: %f \n assym grad: %f\n", proj_grad, fd_grad);
            }
        }
    }
}

Vector3 InverseSolver::find_total_g_grad() {
    Vector3 total_grad = Vector3::zero();
    for (Face f: forwardSolver->hullMesh->faces()){
        total_grad += (goal_area[f] - boundaryBuilder->face_region_area[f]) * 
                        per_face_G_gradient[f];
    }
    return total_grad;
}


// vertex grads
void InverseSolver::find_per_face_per_vertex_grads(bool check_FD) {
    Vector3 zvec = Vector3::zero();
    Vector3 G = forwardSolver->get_G();
    per_face_per_vertex_gradient = FaceData<VertexData<Vector3>>(*forwardSolver->hullMesh);
    for (Face f: forwardSolver->hullMesh->faces()){
        per_face_per_vertex_gradient[f] = VertexData<Vector3>(*forwardSolver->hullMesh, zvec);
        for (Halfedge he: f.adjacentHalfedges()){
            Vertex v0 = he.tailVertex(),
                   v1 = he.tipVertex(),
                   v2 = he.next().tipVertex();
            Vector3 p0 = forwardSolver->hullGeometry->inputVertexPositions[v0], 
                    p1 = forwardSolver->hullGeometry->inputVertexPositions[v1], 
                    p2 = forwardSolver->hullGeometry->inputVertexPositions[v2];
            per_face_per_vertex_gradient[f][v0] += dihedral_angle_grad_B(G, p1, p0, p2);
            per_face_per_vertex_gradient[f][v1] += dihedral_angle_grad_A(G, p1, p0, p2);
            per_face_per_vertex_gradient[f][v2] += dihedral_angle_grad_C(G, p1, p0, p2);
        }
    }
}

VertexData<Vector3> InverseSolver::find_per_vertex_total_grads() {
    Vector3 zvec = Vector3::zero();
    VertexData<Vector3> per_vertex_total_grads(*forwardSolver->hullMesh);
    for (Face f: forwardSolver->hullMesh->faces()){
        for (Vertex v: f.adjacentVertices()){
            per_vertex_total_grads[v] += per_face_per_vertex_gradient[f][v] * 
                                         (goal_area[f] - boundaryBuilder->face_region_area[f]);
        }
    }
    return per_vertex_total_grads;
}

// Vertex grads; Uni mass
void InverseSolver::find_per_face_per_vertex_uni_mass_grads(){
    per_face_per_vertex_uni_mass_gradient = FaceData<VertexData<Vector3>>(*forwardSolver->hullMesh);
    // find
}
        