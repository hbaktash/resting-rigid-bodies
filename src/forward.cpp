#include "forward.h"


ForwardSolver::ForwardSolver(ManifoldSurfaceMesh* inputMesh_, VertexPositionGeometry* inputGeo_){
    mesh = inputMesh_;
    geometry  = inputGeo_;
}


void ForwardSolver::center_geometry(){
    Vector3 center{0.,0.,0.};
    for(Vertex v: mesh->vertices()){
        center += geometry->inputVertexPositions[v];
    }
    center /= (double) mesh->nVertices(); 
    for(Vertex v: mesh->vertices()){
        geometry->inputVertexPositions[v] -= center;
    }
}

void ForwardSolver::align_geometry(Vector2 dir){
    center_geometry();
    double theta = atan2(dir.y, dir.x);
    double to_rotate = 1.5*PI - theta;
    for (Vertex v: mesh->vertices()){
        Vector3 pos = geometry->inputVertexPositions[v];
        Vector2 tmp_v2{pos.x, pos.y};
        tmp_v2 = tmp_v2.rotate(to_rotate);
        geometry->inputVertexPositions[v] = {tmp_v2.x, tmp_v2.y, 0.};
    }
}

void ForwardSolver::find_contact(Vector2 initial_ori){
    align_geometry(initial_ori);
    Vertex contact_point = mesh->vertex(0);
    // choose lowest y
    for (Vertex v: mesh->vertices()){
        Vector3 pos = geometry->inputVertexPositions[v];
        if (pos.y <= geometry->inputVertexPositions[contact_point].y)
            contact_point = v;
    }
    // update state
    Vertex v = contact_point;
    curr_state.first = v;
    curr_state.second = v;
}