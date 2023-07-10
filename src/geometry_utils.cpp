#include "geometry_utils.h"

double signed_volume(Vector3 a, Vector3 b, Vector3 c, Vector3 d){
    return (1.0/6.0)*dot(cross(b-a,c-a),d-a);
}

Vector3 find_center_of_mass(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geometry){
    Vector3 O({0.,0.,0.}); // could be anywhere
    Vector3 G({0.,0.,0.});
    double total_volume = 0.;
    for (Face f: mesh.faces()){
        Halfedge first_he = f.halfedge(),
                 next_he = f.halfedge().next();
        Vertex v0 = first_he.tailVertex();
        Vector3 p0 = geometry.inputVertexPositions[v0];
        while (next_he != first_he){
            Vertex va = next_he.tailVertex(), 
                   vb = next_he.tipVertex();
            Vector3 pa = geometry.inputVertexPositions[va],
                    pb = geometry.inputVertexPositions[vb];
            double curr_signed_vol = signed_volume(O, p0, pa, pb);
            Vector3 tmp_center = (O + p0 + pa + pb)/4.;
            G += tmp_center * curr_signed_vol;
            total_volume += curr_signed_vol;
            next_he = next_he.next();
        }
    }
    return G/total_volume;
}


void center_and_normalize(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry){
    Vector3 G({0.,0.,0.});
    for (Vertex v: mesh->vertices()){
        G += geometry->inputVertexPositions[v];
    }
    G /= (double)mesh->nVertices();
    double max_radi = 0.0;
    for (Vertex v: mesh->vertices()){
        geometry->inputVertexPositions[v] -= G;
        double tmp_norm = norm(geometry->inputVertexPositions[v]);
        if (tmp_norm >= max_radi){
            max_radi = tmp_norm;
        }
    }
    for (Vertex v: mesh->vertices()){
        geometry->inputVertexPositions[v] /= max_radi;
    }
}