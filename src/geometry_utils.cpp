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

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"





// tmp tools ; should be in utils
Vector3 vec2vec3(Vector<double> v){
    Vector3 ans({v[0], v[1], v[2]});
    return ans;
}


Vector3 vec3d_to_vec3(Eigen::Vector3d v){
    Vector3 ans({v[0], v[1], v[2]});
    return ans;
}

// tiny geometry stuff; move elsewhere
// Gradients
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
    return -dihedral_angle_grad_G(A, G, B, C);
}
// wing side vertex B; kinda simple
Vector3 dihedral_angle_grad_B(Vector3 G, Vector3 A, Vector3 B, Vector3 C){
    Vector3 nB_un = cross(A - G, B - A);
    return  nB_un * (G-A).norm()/nB_un.norm2();
}
// wing side vertex C; kinda simple
Vector3 dihedral_angle_grad_C(Vector3 G, Vector3 A, Vector3 B, Vector3 C){
    Vector3 nC_un = cross(A - G, A - C);
    return nC_un * (G-A).norm()/nC_un.norm2();
}


// big geometry stuff


double signed_volume(Vector3 a, Vector3 b, Vector3 c, Vector3 d){
    return (1.0/6.0)*dot(cross(b-a,c-a),d-a);
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

double polygonal_face_area(Face f, VertexPositionGeometry &geometry){
    Halfedge first_he = f.halfedge();
    Vertex v0 = first_he.tailVertex(),
           v1 = first_he.tipVertex();
    Vector3 p0 = geometry.inputVertexPositions[v0],
            p1 = geometry.inputVertexPositions[v1];
    double area = 0.;
    for (Halfedge he: f.adjacentHalfedges()){
        if(he != first_he && he.tipVertex() != v0){
            Vertex v2 = he.tipVertex();
            Vector3 p2 = geometry.inputVertexPositions[v2];
            area += 0.5 * norm(cross(p1 - p0, p2 - p0));
        }
    }
    return area;
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
    if (total_volume < 0){
        throw std::logic_error("total vol < 0; proly bad normal orientation\n");
        // total_volume *= -1.;
        // for (Face f: mesh.faces())
        //     mesh.invertOrientation(f);
    }

    return {G/total_volume, total_volume};
}


void center_and_normalize(SurfaceMesh* mesh, VertexPositionGeometry* geometry){
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

bool check_hollow_tet_vertex(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry, Vertex v){
    printf("      at v %d deg: %d\n", v.getIndex(), v.degree());
    if (v.degree() > 3) return false;
    if (v.degree() < 3) throw std::logic_error("what the fuck?\n");
    size_t count = 0;
    for (Edge e: v.adjacentEdges()){
        double dihedangle = geometry->edgeDihedralAngle(e);
        printf("      at v %d diherdal angle: %f\n", v.getIndex(), dihedangle);
        if (dihedangle < 0. || dihedangle > PI - 1e-4)
            count++;
    }
    return count == 3;
}

bool check_convexity_and_repair(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry){
    // mesh->removeVertex()
    // NOTE: has to be triangle mesh, Manifold, no degree 1 vertex
    size_t e_count = 0, trials = 0;
    printf("$$$$$ here\n");
    while (e_count < mesh->nEdges() - 1){
        e_count = 0;
        for (Edge e: mesh->edges()){ // TODO: how to handle this?
            double dihedangle = geometry->edgeDihedralAngle(e);
            // double normals_dot = dot(geometry->faceNormal(e.halfedge().face()), // assuming outward normals
                                    //  geometry->faceNormal(e.halfedge().twin().face()));
            if (dihedangle < 0.){
                printf("$$$$$ here for edge %d: %d, %d angle %f, $$$$\n", e.getIndex(), e.firstVertex().getIndex(), e.secondVertex().getIndex(), dihedangle);
                if (check_hollow_tet_vertex(mesh, geometry, e.firstVertex())){
                    printf("got a vertex! %d\n", e.firstVertex().getIndex());
                    Face new_face = mesh->removeVertex(e.firstVertex());
                    mesh->compress();
                    break; // since edges get removed
                }
                else if (check_hollow_tet_vertex(mesh, geometry, e.secondVertex())){
                    printf("got a vertex! %d\n", e.secondVertex().getIndex());
                    Face new_face = mesh->removeVertex(e.secondVertex());
                    mesh->compress();
                    break; // since edges get removed
                }
                else {
                    mesh->flip(e);
                    mesh->compress();
                    break;
                }
            }
            else e_count++;
        }
        // mesh->compress();
        trials++;
    }
    return trials > 1;
}


Edge single_convexity_repair(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry){
    // mesh->removeVertex()
    // NOTE: has to be triangle mesh, Manifold, no degree 1 vertex
    size_t e_count = 0;
    printf("$$$$$ here\n");
    // while (e_count < mesh->nEdges() - 1){
    e_count = 0;
    for (Edge e: mesh->edges()){ // TODO: how to handle this?
        double dihedangle = geometry->edgeDihedralAngle(e);
        // double normals_dot = dot(geometry->faceNormal(e.halfedge().face()), // assuming outward normals
                                //  geometry->faceNormal(e.halfedge().twin().face()));
        if (dihedangle < 0. || dihedangle > PI - 1e-4){
            printf("$$$$$ here for edge %d: %d, %d angle %f, $$$$\n", e.getIndex(), e.firstVertex().getIndex(), e.secondVertex().getIndex(), dihedangle);
            printf("$$$$$ diherdal angle: %f\n", dihedangle);
            if (check_hollow_tet_vertex(mesh, geometry, e.firstVertex())){
                printf("got a vertex! %d\n", e.firstVertex().getIndex());
                Face new_face = mesh->removeVertex(e.firstVertex());
                mesh->compress();
                // break; // since edges get removed
            }
            else if (check_hollow_tet_vertex(mesh, geometry, e.secondVertex())){
                printf("got a vertex! %d\n", e.secondVertex().getIndex());
                Face new_face = mesh->removeVertex(e.secondVertex());
                mesh->compress();
                // break; // since edges get removed
            }
            else {
                if (e.getIndex() != INVALID_IND){
                    Vertex v1 = e.firstVertex(),
                            v2 = e.secondVertex(),
                            v3 = e.halfedge().next().tipVertex(),
                            v4 = e.halfedge().twin().next().tipVertex();
                    std::vector<std::vector<size_t>> ff({{0, 1, 2}, {0,1,3}});
                    std::vector<Vector3> possss({geometry->inputVertexPositions[v1],
                                                 geometry->inputVertexPositions[v2],
                                                 geometry->inputVertexPositions[v3],
                                                 geometry->inputVertexPositions[v4]});
                    polyscope::registerSurfaceMesh("temp old mesh", possss,
                                             ff);  
                }
                mesh->flip(e);
                mesh->compress();
                //DEBUG
                return e;
                // break;
            }
            return Edge();
        }
        else e_count++;
    }
        // mesh->compress();
    // }
}


DenseMatrix<double> solve_dense_b(LinearSolver<double> *solver, DenseMatrix<double> b){
    DenseMatrix<double> sol(b.rows(), b.cols());
    
    for (size_t i = 0; i < b.cols(); i++){
        sol.col(i) = solver->solve(b.col(i));
    }

    return sol;
}


Eigen::MatrixXd sobolev_diffuse_gradients(Eigen::MatrixXd grads, ManifoldSurfaceMesh &hull_mesh, 
                                              double sobolev_lambda, size_t sobolev_p){
    if (sobolev_lambda == 0.){
        return grads;
    }
    size_t n = hull_mesh.nVertices(), 
           e = hull_mesh.nEdges();
    SparseMatrix<double> graph_L(n, n);
    std::vector<Eigen::Triplet<double>> gL_tripletList;
    gL_tripletList.reserve(4*e);
    for (Edge e: hull_mesh.edges()){
        Vertex v1 = e.halfedge().vertex(),
               v2 = e.halfedge().twin().vertex();
        gL_tripletList.push_back(Eigen::Triplet<double>(v1.getIndex(), v2.getIndex(), -1.));
        gL_tripletList.push_back(Eigen::Triplet<double>(v2.getIndex(), v1.getIndex(), -1.));
        gL_tripletList.push_back(Eigen::Triplet<double>(v1.getIndex(), v1.getIndex(), 1.));
        gL_tripletList.push_back(Eigen::Triplet<double>(v2.getIndex(), v2.getIndex(), 1.));
    }
    graph_L.setFromTriplets(gL_tripletList.begin(), gL_tripletList.end());

    // Sobolev operator
    SparseMatrix<double> sobolevOp = sobolev_lambda * graph_L + identityMatrix<double>(n);
    PositiveDefiniteSolver<double> sobolevSolver(sobolevOp);
    DenseMatrix<double> sobolev_grads;
    for (size_t i = 0; i < sobolev_p; i++){
        sobolev_grads = solve_dense_b(&sobolevSolver, grads);
    }
    
    return sobolev_grads;                                                    
}


VertexData<Vector3> sobolev_diffuse_gradients(VertexData<Vector3> grads, ManifoldSurfaceMesh &hull_mesh, double sobolev_lambda, size_t sobolev_p){
    Eigen::MatrixXd sobolev_grads = sobolev_diffuse_gradients(vertex_data_to_matrix(grads), hull_mesh, sobolev_lambda, sobolev_p);
    VertexData<Vector3> sobolev_grads_vd(hull_mesh);
    for (Vertex v: hull_mesh.vertices()){
        sobolev_grads_vd[v] = vec2vec3(sobolev_grads.row(v.getIndex()));
    }
    return sobolev_grads_vd;
}