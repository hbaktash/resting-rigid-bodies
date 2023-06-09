#include "forward3D.h"


Forward3DSolver::Forward3DSolver(ManifoldSurfaceMesh* inputMesh_, VertexPositionGeometry* inputGeo_,
                             Vector3 inputG_){
    mesh = inputMesh_;
    geometry = inputGeo_;
    G = inputG_;
}


// just the Gaussian curvature
void Forward3DSolver::compute_vertex_probabilities(){
    vertex_probabilities = VertexData<double>(*hullMesh, 0.);
    for (Vertex v: hullMesh->vertices()){
        vertex_probabilities[v] = hullGeometry->vertexGaussianCurvature(v)/(4*PI);
    }
}

// stabilizable := the normal from G can touch the ground
// stable := the normal from G falls withing the element
bool Forward3DSolver::vertex_is_stablizable(Vertex v){
    Vector3 Gp = hullGeometry->inputVertexPositions[v] - G;
    double gp_norm = dot(Gp, Gp);
    for (Vertex other_v: hullMesh->vertices()){
        if (v != other_v){
            Vector3 Gp2 = hullGeometry->inputVertexPositions[other_v] - G;
            if (dot(Gp2, Gp) > gp_norm)
                return false;
        }
    }
    return true;
}
bool Forward3DSolver::edge_is_stable(Edge e){
    Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
    return false;
}
bool Forward3DSolver::edge_is_stablizable(Edge e){
    return false;
}
bool Forward3DSolver::face_is_stable(Face f){
    return false;
}