#include "markov_model.h"



// trivial constructors
SudoFace::SudoFace(Edge host_edge_, Vector3 normal_, SudoFace *next_sudo_face_, SudoFace *prev_sudo_face_){
    host_edge = host_edge_;
    normal = normal_;
    next_sudo_face = next_sudo_face_;
    prev_sudo_face = prev_sudo_face_;
}

// SudoEdge::SudoEdge(Edge host_edge_, SudoFace *first_sudo_face_, SudoFace *second_sudo_face_){
//     host_edge = host_edge_;
//     first_sudo_face = first_sudo_face_;
//     second_sudo_face = second_sudo_face_;
// }

RollingMarkovModel::RollingMarkovModel(Forward3DSolver *forward_solver_){
    forward_solver = forward_solver_;
    G = forward_solver->G;
    mesh = forward_solver->hullMesh;
    geometry = forward_solver->hullGeometry;
}

RollingMarkovModel::RollingMarkovModel(ManifoldSurfaceMesh* mesh_, VertexPositionGeometry* geometry_, Vector3 G_){
    G = G_;
    mesh = mesh_;
    geometry = geometry_;
    forward_solver = new Forward3DSolver(mesh, geometry, G);
}


// split the SudoEdge starting with this SudoFace 
// TODO: Assign source/sink while/before/after finding the new_normal
SudoFace* SudoFace::split_sudo_edge(Vector3 new_normal){
    // current SudoEdge will be the first, by contract
    if (this == next_sudo_face){
        printf("This is a terminal SudoFace");
        return nullptr;
    }
    // SudoEdge is well-defined 
    SudoFace *new_sudo_face = new SudoFace(host_edge, new_normal, next_sudo_face, this);
    next_sudo_face->prev_sudo_face = new_sudo_face;
    this->next_sudo_face = new_sudo_face;
    return new_sudo_face;
}