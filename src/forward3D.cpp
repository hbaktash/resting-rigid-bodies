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
    Vector3 p1 = hullGeometry->inputVertexPositions[v1], p2 = hullGeometry->inputVertexPositions[v2];
    if (dot(G - p1, p2-p1) <= 0 ||
        dot(G - p2, p1-p2) <= 0) // Gv1v2 triange is obtuse
        return false;
    return true;
}


// Vector3 


bool Forward3DSolver::edge_is_stablizable(Edge e){
    // finding the orthogonal vector from G onto AB 
    Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
    Vector3 A = hullGeometry->inputVertexPositions[v1], B = hullGeometry->inputVertexPositions[v2];
    Vector3 GB = B - G,
            AB = B - A;
    Vector3 ortho_g = GB - AB*dot(AB, GB)/dot(AB,AB);
    tmp_test_vec = ortho_g;
    // checking if other vertices are more aligned with orhto_G
    double ortho_g_norm = dot(ortho_g, ortho_g);
    for (Vertex other_v: hullMesh->vertices()){
        if (other_v != e.firstVertex() && other_v != e.secondVertex()){
            Vector3 GP = hullGeometry->inputVertexPositions[other_v] - G;
            if (dot(GP,ortho_g) > ortho_g_norm){
                return false;
            }
        }
    }
    return true;
}


bool Forward3DSolver::face_is_stable(Face f){
    // need to check the 3 dihedral angles
    // build the base triangle
    Halfedge he = f.halfedge();
    Vertex v1 = he.tailVertex(), v2 = he.tipVertex(), v3 = he.next().tipVertex();
    assert(v1 == he.next().next().tipVertex());
    Vector3 A = hullGeometry->inputVertexPositions[v1], 
            B = hullGeometry->inputVertexPositions[v2], 
            C = hullGeometry->inputVertexPositions[v3];
    // build compatible normals; dir: A->B->C
    Vector3 N_base = cross(B - A, C - B),
            N_GAB  = cross(A - B, G - A),
            N_GBC  = cross(B - C, G - B),
            N_GCA  = cross(C - A, G - C);
    if (dot(N_base, N_GAB) >= 0 ||
        dot(N_base, N_GBC) >= 0 ||
        dot(N_base, N_GCA) >= 0)
        return false;
    return true;
}