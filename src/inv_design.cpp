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

InverseSolver::InverseSolver(BoundaryBuilder* boundaryBuilder){
    this->boundaryBuilder = boundaryBuilder;
    forwardSolver = boundaryBuilder->forward_solver;
    set_fair_distribution();
}

void InverseSolver::set_fair_distribution() {
    goal_area = FaceData<double>(*forwardSolver->hullMesh, 
                                         1./(double)forwardSolver->hullMesh->nFaces());
}


void InverseSolver::find_per_face_G_grads(){
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
}

Vector3 InverseSolver::find_total_g_grad() {
    Vector3 total_grad = Vector3::zero();
    for (Face f: forwardSolver->hullMesh->faces()){
        total_grad += (goal_area[f] - boundaryBuilder->face_region_area[f]) * 
                        per_face_G_gradient[f];
    }
    return total_grad;
}