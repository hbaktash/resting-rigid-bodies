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

#include "geometry_utils.h"

double signed_volume(Vector3 a, Vector3 b, Vector3 c, Vector3 d){
    return (1.0/6.0)*dot(cross(b-a,c-a),d-a);
}

// check if G is inside the polyehdra (positive mass)
bool G_is_inside(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geometry, Vector3 G){
    return true;
    for (Face f : mesh.faces()){
        // assuming planar faces
        Vertex v = f.halfedge().vertex();
        Vector3 p = geometry.inputVertexPositions[v];
        // assuming outward face normals
        if(dot(p - G, geometry.faceNormal(f)) <= 0.)
            return false;
    }
    return true;
}

std::pair<Vector3, double> find_center_of_mass(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geometry){
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
    if (total_volume < 0)
        throw std::logic_error("total vol < 0; proly bad normal orientation\n");
    return {G/total_volume, total_volume};
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

double area(Vector3 A, Vector3 B, Vector3 C){
    return 0.5*norm(cross(B-A, C-A));
}

double ray_intersect_triangle(Vector3 O, Vector3 v, Vector3 A, Vector3 B, Vector3 C){
    Vector3 e1 = B - A,
            e2 = C - A;
    Vector3 n = cross(e1, e2);
    if (dot(v,n) < 1e-8) return -1;
    double t = (dot(A - O,n))/dot(v,n);
    if (t < 0)
        return -1;

    Vector3 P = O + t*v;
    if (area(A,B,C) >= area(P,A,B) + area(P,A,C)+ area(P,C,B) - 1e-6)
        return t;
    else 
        return -1;
}

double ray_intersect(Vector3 O, Vector3 v, std::vector<Vector3> polygon){
    size_t n = polygon.size();
    for (size_t i = 1; i < n; i++){
        double tmp_t = ray_intersect_triangle(O, v, polygon[0], polygon[i], polygon[(i+1)%n]);
        if (tmp_t != -1)
            return tmp_t;
    }
    return -1;
}