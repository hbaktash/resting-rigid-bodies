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

// optimization stuff

InverseSolver::InverseSolver(BoundaryBuilder* boundaryBuilder){
    this->boundaryBuilder = boundaryBuilder;
    forwardSolver = boundaryBuilder->forward_solver;
    set_fair_distribution();
}

void InverseSolver::set_fair_distribution() {
    goal_prob = FaceData<double>(*forwardSolver->hullMesh, 
                                         1./(double)forwardSolver->hullMesh->nFaces());
}

void InverseSolver::set_fair_distribution_for_sink_faces(){
    size_t count = 0;
    goal_prob = FaceData<double>(*forwardSolver->hullMesh, 0.);
    for (Face f: forwardSolver->hullMesh->faces()){
        if (flow_structure[f] == f) // f is sink
            count++;
    }
    for (Face f: forwardSolver->hullMesh->faces()){
        if (flow_structure[f] == f){
            goal_prob[f] = 1./(double)count;
            old_stable_normals.push_back(forwardSolver->hullGeometry->faceNormal(f));
        }
    }
}

void InverseSolver::update_fair_distribution(double normal_threshold){
    printf("updating fair dist\n");
    goal_prob = FaceData<double>(*forwardSolver->hullMesh, 0.);
    size_t old_count = old_stable_normals.size();
    std::vector<Face> new_stable_faces, all_stable_faces;
    for (Face f: forwardSolver->hullMesh->faces()){
        if (flow_structure[f] == f){
            all_stable_faces.push_back(f);
        }
    }
    std::vector<Vector3> new_stable_normals;
    size_t new_stables_count = 0;
    for(Vector3 old_stable_normal: old_stable_normals){
        double closest_normal_dist = 10.; // anything < 2 would work
        Face closest_face = Face();
        Vector3 closest_f_normal;
        for (Face f: all_stable_faces){
            Vector3 f_normal = forwardSolver->hullGeometry->faceNormal(f);
            if ((old_stable_normal - f_normal).norm() < closest_normal_dist){
                closest_face = f;
                closest_f_normal = f_normal;
                closest_normal_dist = (old_stable_normal - f_normal).norm();
            }
        }
        if (closest_normal_dist < normal_threshold){
            printf("    - updated to face %d\n", closest_face.getIndex());
            new_stable_faces.push_back(closest_face);
            new_stables_count++;
            new_stable_normals.push_back(closest_f_normal);
        }
        else 
            printf("    - face %d removed due to threshold %f\n", closest_face.getIndex(), normal_threshold);
    }
    double goal_fair_prob = 1./(double)new_stables_count;
    for (Face f: new_stable_faces)
        goal_prob[f] = goal_fair_prob;
    old_stable_normals = new_stable_normals;
    printf("done updating stable faces\n");
}

void InverseSolver::find_d_pf_d_Gs(bool check_FD){
    Vector3 zero_vec = Vector3::zero();
    Vector3 G = forwardSolver->get_G();
    d_pf_d_G = FaceData<Vector3>(*forwardSolver->hullMesh, zero_vec);
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
        d_pf_d_G[f] = f_g_grad;
    }
    if (check_FD){
        printf("FD check: \n");
        double step = 1e-6;
        Vector3 e_x({1., 0. ,0.}),
                e_y({0., 1. ,0.}),
                e_z({0., 0. ,1.});
        boundaryBuilder->build_boundary_normals();
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
                double proj_grad = dot(d_pf_d_G[f], dG),
                       fd_grad = (new_face_areas[f.getIndex()] - old_face_areas[f.getIndex()])/step;
                if (proj_grad - fd_grad > 1e-4){
                    printf("Missmatch:\n");
                    printf(" proj grad: %f \n assym grad: %f\n", proj_grad, fd_grad);
                }
            }
        }
        forwardSolver->set_G(G);
        forwardSolver->initialize_pre_computes();   
    }
}

Vector3 InverseSolver::find_total_g_grad() {
    Vector3 total_grad = Vector3::zero();
    for (Face f: forwardSolver->hullMesh->faces()){
        total_grad += (goal_prob[f] * 4.*PI - boundaryBuilder->face_region_area[f]) * 
                        d_pf_d_G[f];
    }
    return total_grad;
}


// vertex grads
void InverseSolver::find_d_pf_dvs(bool check_FD) {
    if (!forwardSolver->updated)
        forwardSolver->initialize_pre_computes();
    Vector3 zvec = Vector3::zero();
    Vector3 G = forwardSolver->get_G();
    d_pf_dv = FaceData<VertexData<Vector3>>(*forwardSolver->hullMesh);
    for (Face f: forwardSolver->hullMesh->faces()){
        d_pf_dv[f] = VertexData<Vector3>(*forwardSolver->hullMesh, zvec);
        for (Halfedge he: f.adjacentHalfedges()){
            Vertex v0 = he.tailVertex(),
                   v1 = he.tipVertex(),
                   v2 = he.next().tipVertex();
            Vector3 p0 = forwardSolver->hullGeometry->inputVertexPositions[v0], 
                    p1 = forwardSolver->hullGeometry->inputVertexPositions[v1], 
                    p2 = forwardSolver->hullGeometry->inputVertexPositions[v2];
            d_pf_dv[f][v0] += dihedral_angle_grad_B(G, p1, p0, p2);
            d_pf_dv[f][v1] += dihedral_angle_grad_A(G, p1, p0, p2);
            d_pf_dv[f][v2] += dihedral_angle_grad_C(G, p1, p0, p2);
        }
    }
    if (check_FD){
        printf("FD check: \n");
        double step = 1e-6;
        boundaryBuilder->build_boundary_normals();
        Vector<double> old_face_areas = boundaryBuilder->face_region_area.toVector();
        for (Vertex v: forwardSolver->hullMesh->vertices()){
            printf("- at v %d\n", v.getIndex());
            Vector3 old_p = forwardSolver->hullGeometry->inputVertexPositions[v];
            Vector3 e_x({1., 0. ,0.}),
                    e_y({0., 1. ,0.}),
                    e_z({0., 0. ,1.});
            for (Vector3 dp: {e_x, e_y, e_z}){
                Vector3 tmp_p = old_p + dp * step;
                forwardSolver->hullGeometry->inputVertexPositions[v] = tmp_p;
                forwardSolver->updated = false;
                forwardSolver->initialize_pre_computes();
                boundaryBuilder->build_boundary_normals();
                Vector<double> new_face_areas = boundaryBuilder->face_region_area.toVector();
                for (Face f: v.adjacentFaces()){
                    double proj_grad = dot(d_pf_dv[f][v], dp),
                           fd_grad = (new_face_areas[f.getIndex()] - old_face_areas[f.getIndex()])/step;
                    if (proj_grad - fd_grad > 1e-4){
                        printf("Missmatch:\n");
                        std::cout << " - at dp: "<< dp <<"\n";
                        printf("   - proj grad: %f \n   - assym grad: %f\n", proj_grad, fd_grad);
                    }
                }
                // break;
            }
            forwardSolver->hullGeometry->inputVertexPositions[v] = old_p;
            forwardSolver->updated = false;
            forwardSolver->initialize_pre_computes();   
            // break; 
        }
    }
}

VertexData<Vector3> InverseSolver::find_total_vertex_grads() {
    Vector3 zvec = Vector3::zero();
    VertexData<Vector3> per_vertex_total_grads(*forwardSolver->hullMesh, zvec);
    for (Face f: forwardSolver->hullMesh->faces()){
        printf("   at f %d, prob: %f goal area: %f \n", f.getIndex(), boundaryBuilder->face_region_area[f]/(4.*PI),
                                                                      goal_prob[f]);
        for (Vertex v: f.adjacentVertices()){
            per_vertex_total_grads[v] += d_pf_dv[f][v] * 
                                         (goal_prob[f] * 4 * PI - boundaryBuilder->face_region_area[f]);
        }
    }
    return per_vertex_total_grads;
}


Vector<double> vec32vec(Vector3 v){
    Vector<double> ans(3);
    ans[0] = v.x;
    ans[1] = v.y;
    ans[2] = v.z;
    return ans;
}

// populating VertexData<DenseMatrix<double>> dv_d_G;
void InverseSolver::find_dG_dvs(){
    // TODO; generalize for polygonal faces
    if (!forwardSolver->updated)
        forwardSolver->initialize_pre_computes();
    Vector3 G = forwardSolver->get_G();
    DenseMatrix<double> zmat = DenseMatrix<double>::Zero(3,3);
    dG_dv = VertexData<DenseMatrix<double>>(*forwardSolver->hullMesh, zmat);
    // for (Vertex v: forwardSolver->hullMesh->vertices())
    //     dG_dv[v] = zmat;
    for (Face f: forwardSolver->hullMesh->faces()){
        // double face_area = forwardSolver->inputGeometry->faceArea(f);
        double face_area = polygonal_face_area(f, *forwardSolver->hullGeometry);

        Vector3 face_normal = forwardSolver->inputGeometry->faceNormal(f); // assuming outward normals
        size_t face_degree = f.degree();
        // assuming polygon faces here; TODO; check things
        Vector3 vert_sum = Vector3::zero();
        for (Vertex tmp_v: f.adjacentVertices())
            vert_sum += forwardSolver->inputGeometry->inputVertexPositions[tmp_v];
        for (Halfedge he: f.adjacentHalfedges()){
            Vertex v = he.tailVertex();
            Vector3 p = forwardSolver->hullGeometry->inputVertexPositions[v];
            Vector3 Gf_G = (vert_sum + p)/(double)(face_degree + 1) - G;
            DenseMatrix<double> tmp_mat = vec32vec(Gf_G) * 
                                          vec32vec(face_normal).transpose();
            assert(tmp_mat.cols() == 3);
            assert(tmp_mat.rows() == 3);
            dG_dv[v] += face_area * tmp_mat;
        }
    }
    double volume = forwardSolver->volume;
    for (Vertex v: forwardSolver->hullMesh->vertices())
        dG_dv[v] /= 3.*volume;
}

// Vertex grads; Uni mass
void InverseSolver::find_uni_mass_d_pf_dv(bool check_FD){
    //  assuming that we only move convex hull vertices
    //  and that input is convex
    find_d_pf_d_Gs();
    find_d_pf_dvs();
    find_dG_dvs();
    uni_mass_d_pf_dv = FaceData<VertexData<Vector3>>(*forwardSolver->hullMesh);
    Vector3 zvec = Vector3::zero();
    Vector3 G = forwardSolver->get_G();
    for (Face f: forwardSolver->hullMesh->faces()){
        uni_mass_d_pf_dv[f] = VertexData<Vector3>(*forwardSolver->hullMesh, zvec);
        for (Vertex v: f.adjacentVertices()){
            // if input not convex, just map hullMesh to input mesh 
            // Vector3 tmp_d_pf_dv = d_pf_dv[f][v];
            Vector<double> term2 = dG_dv[v].transpose() * vec32vec(d_pf_d_G[f]);
            uni_mass_d_pf_dv[f][v] = d_pf_dv[f][v] + Vector3({term2[0], term2[1], term2[2]});
        }
    }
    if (check_FD){
        printf("FD check: \n");
        double step = 1e-6;
        forwardSolver->set_uniform_G();
        forwardSolver->initialize_pre_computes();        
        boundaryBuilder->build_boundary_normals();
        Vector<double> old_face_areas = boundaryBuilder->face_region_area.toVector();

        for (Vertex v: forwardSolver->hullMesh->vertices()){
            printf("- at v %d\n", v.getIndex());
            Vector3 old_p = forwardSolver->hullGeometry->inputVertexPositions[v];
            Vector3 e_x({1., 0. ,0.}),
                    e_y({0., 1. ,0.}),
                    e_z({0., 0. ,1.});
            for (Vector3 dp: {e_x, e_y, e_z}){
                Vector3 tmp_p = old_p + dp * step;
                forwardSolver->hullGeometry->inputVertexPositions[v] = tmp_p;
                forwardSolver->set_uniform_G();
                forwardSolver->initialize_pre_computes();
                boundaryBuilder->build_boundary_normals();
                Vector<double> new_face_areas = boundaryBuilder->face_region_area.toVector();
                for (Face f: v.adjacentFaces()){
                    double proj_grad = dot(uni_mass_d_pf_dv[f][v], dp),
                           fd_grad = (new_face_areas[f.getIndex()] - old_face_areas[f.getIndex()])/step;
                    if (proj_grad - fd_grad > 1e-4){
                        printf("Missmatch:\n");
                        std::cout << " - at dp: "<< dp <<"\n";
                        printf("   - proj grad: %f \n assym grad: %f\n", proj_grad, fd_grad);
                    }
                }
                // break;
            }
            forwardSolver->hullGeometry->inputVertexPositions[v] = old_p;
            forwardSolver->updated = false;
            forwardSolver->initialize_pre_computes();   
            // break; 
        }
    }
}

// 
VertexData<Vector3> InverseSolver::find_uni_mass_total_vertex_grads(bool with_flow_structure,
                                                                    double stable_normal_update_thres){
    Vector3 zvec = Vector3::zero();
    VertexData<Vector3> uni_mass_total_vertex_grads(*forwardSolver->hullMesh, zvec);
    if (with_flow_structure && stable_normal_update_thres < 0)
        set_fair_distribution_for_sink_faces();
    else if (stable_normal_update_thres > 0)
        update_fair_distribution(stable_normal_update_thres);
    for (Face f: forwardSolver->hullMesh->faces()){
        Face last_sink_face;
        for (Vertex v: f.adjacentVertices()){
            if (with_flow_structure) {
                last_sink_face = flow_structure[f];
                uni_mass_total_vertex_grads[v] += uni_mass_d_pf_dv[f][v] * 
                                            (goal_prob[last_sink_face]* 4.*PI - boundaryBuilder->face_region_area[last_sink_face]);
            }
            else {
                uni_mass_total_vertex_grads[v] += uni_mass_d_pf_dv[f][v] * 
                                            (goal_prob[f]* 4.*PI - boundaryBuilder->face_region_area[f]);
            }
        }
        printf("   finding total grad for f %d area: %f goal area: %f \n", f.getIndex(),
                                                                        boundaryBuilder->face_region_area[last_sink_face]/(4.*PI),
                                                                        goal_prob[last_sink_face]);
    }
    return uni_mass_total_vertex_grads;
}