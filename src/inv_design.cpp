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
    std::vector<double> face_areas;
    size_t goal_stable_count = 6; // TODO: take as input?
    size_t count = goal_stable_count;
    goal_prob = FaceData<double>(*forwardSolver->hullMesh, 0.);
    for (Face f: forwardSolver->hullMesh->faces()){
        if (flow_structure[f] == f) { // f is sink
            face_areas.push_back(boundaryBuilder->face_region_area[f]);
        }
    }
    // select the high prob faces for optimizations
    std::sort(face_areas.begin(), face_areas.end());
    double nth_largest_area = face_areas[face_areas.size() - goal_stable_count];
    for (Face f: forwardSolver->hullMesh->faces()){
        Vector3 face_normal = forwardSolver->hullGeometry->faceNormal(f);
        if (flow_structure[f] == f && boundaryBuilder->face_region_area[f] >= nth_largest_area){
            goal_prob[f] = 1./(double)count;
            old_stable_normals.push_back(face_normal);
        }
    }
}


//
void InverseSolver::initialize_interior_vertex_trackers(){
    for (Face f: forwardSolver->hullMesh->faces()){
        Vector3 face_normal = forwardSolver->hullGeometry->faceNormal(f);
        old_normals.push_back(face_normal);
    }
    
    interior_v_to_hull_f = VertexData<Face>(*forwardSolver->inputMesh, Face());
    interior_v_to_hull_f_hit_ratio = VertexData<double>(*forwardSolver->inputMesh, -1.);

    Vector3 O = forwardSolver->get_G();
    for (Vertex v: forwardSolver->inputMesh->vertices()){
        if (forwardSolver->on_hull_index[v] != INVALID_IND)
            continue;
        Vector3 p = forwardSolver->inputGeometry->inputVertexPositions[v];
        for (Face hull_f: forwardSolver->hullMesh->faces()){
            // assuming planar faces
            Halfedge he = hull_f.halfedge();
            std::vector<Vector3> polygon_points;
            for (Vertex fv: hull_f.adjacentVertices()){
                Vector3 p1 = forwardSolver->hullGeometry->inputVertexPositions[fv];
                polygon_points.push_back(p1);
            }
            double t = ray_intersect(O, (p - O).normalize(), polygon_points);
            if (t != -1){
                interior_v_to_hull_f[v] = hull_f;
                interior_v_to_hull_f_hit_ratio[v] = (p - O).norm()/t;
            }
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
    // match with new normals
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

Vector3 vec2vec3(Vector<double> v){
    Vector3 ans({v[0], v[1], v[2]});
    return ans;
}

// populating VertexData<DenseMatrix<double>> dv_d_G;
void InverseSolver::find_dG_dvs(){
    // TODO; generalize for polygonal faces
    if (!forwardSolver->updated)
        forwardSolver->initialize_pre_computes();
    Vector3 G = forwardSolver->get_G();
    DenseMatrix<double> zmat = DenseMatrix<double>::Zero(3,3);
    dG_dv = VertexData<DenseMatrix<double>>(*forwardSolver->inputMesh, zmat);
    for (Face f: forwardSolver->inputMesh->faces()){
        // double face_area = forwardSolver->inputGeometry->faceArea(f);
        double face_area = polygonal_face_area(f, *forwardSolver->inputGeometry);

        Vector3 face_normal = forwardSolver->inputGeometry->faceNormal(f); // assuming outward normals
        size_t face_degree = f.degree();
        // assuming polygon faces here; TODO; check things
        Vector3 vert_sum = Vector3::zero();
        for (Vertex tmp_v: f.adjacentVertices())
            vert_sum += forwardSolver->inputGeometry->inputVertexPositions[tmp_v];
        for (Halfedge he: f.adjacentHalfedges()){
            Vertex v = he.tailVertex();
            Vector3 p = forwardSolver->inputGeometry->inputVertexPositions[v];
            Vector3 Gf_G = (vert_sum + p)/(double)(face_degree + 1) - G;
            DenseMatrix<double> tmp_mat = vec32vec(Gf_G) * 
                                          vec32vec(face_normal).transpose();
            assert(tmp_mat.cols() == 3);
            assert(tmp_mat.rows() == 3);
            dG_dv[v] += face_area * tmp_mat;
        }
    }
    double volume = forwardSolver->volume;
    for (Vertex v: forwardSolver->inputMesh->vertices())
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
            // need to find the index in the original mesh (not hull)
            Vertex org_vertex = forwardSolver->inputMesh->vertex(forwardSolver->org_hull_indices[v]); // converting to vertex for assertion
            // dG_dv is fetched from original mesh, not the convex hull; everything else is fetched from hull faces and geometry
            Vector<double> term2 = dG_dv[org_vertex].transpose() * vec32vec(d_pf_d_G[f]);
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
    if (with_flow_structure && stable_normal_update_thres < 0){
        set_fair_distribution_for_sink_faces(); // 
    }
    else if (stable_normal_update_thres > 0)
        update_fair_distribution(stable_normal_update_thres);

    initialize_interior_vertex_trackers(); // follow inner vertices
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
            // printf("      -- tmp grad norm for v %d: %f \n", v.getIndex(), uni_mass_d_pf_dv[f][v].norm());
        }
        printf("   finding total grad for f %d area: %f goal area: %f \n", f.getIndex(),
                                                                        boundaryBuilder->face_region_area[last_sink_face]/(4.*PI),
                                                                        goal_prob[last_sink_face]);
    }
    return uni_mass_total_vertex_grads;
}

DenseMatrix<double> solve_dense_b(LinearSolver<double> *solver, DenseMatrix<double> b){
    DenseMatrix<double> sol(b.rows(), b.cols());
    
    for (size_t i = 0; i < b.cols(); i++){
        sol.col(i) = solver->solve(b.col(i));
    }

    return sol;
}

VertexData<DenseMatrix<double>> 
InverseSolver::find_rotations(DenseMatrix<double> old_pos, DenseMatrix<double> new_pos){
    VertexData<DenseMatrix<double>> rotations(*forwardSolver->inputMesh);
    SparseMatrix<double> L = forwardSolver->inputGeometry->cotanLaplacian; // assuming required already called
    for (Vertex v: forwardSolver->inputMesh->vertices()){
        size_t n_v = v.degree();
        Vector<size_t> neigh_inds(n_v);
        Vector<double> weights(n_v);
        size_t cnt = 0;
        for (Vertex fv: v.adjacentVertices()){
            neigh_inds(cnt) = fv.getIndex();
            weights(cnt) = -L.coeffRef(v.getIndex(), fv.getIndex());
            cnt++;
        }
        DenseMatrix<double> weights_diag = weights.asDiagonal();
        DenseMatrix<double> P_v = old_pos(v.getIndex(), Eigen::all).replicate(n_v, 1) - 
                                   old_pos(neigh_inds, Eigen::all),
                            P_v_new = new_pos(v.getIndex(), Eigen::all).replicate(n_v, 1) - 
                                       new_pos(neigh_inds, Eigen::all);
        DenseMatrix<double> S_v = P_v.transpose() * weights_diag * P_v_new;

        auto bdcSVD = S_v.bdcSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
        DenseMatrix<double> U = bdcSVD.matrixU(),
                            V = bdcSVD.matrixV();
        DenseMatrix<double> R_v = V * U.transpose();
        rotations[v] = R_v;
    }
    return rotations;
}

DenseMatrix<double> vertex_data_to_matrix(VertexData<Vector3> positions){
    size_t n = positions.size();
    DenseMatrix<double> mat(n, 3);
    for (Vertex v: positions.getMesh()->vertices()){
        Vector3 p = positions[v];
        mat.row(v.getIndex()) = vec32vec(p);
    }
    return mat;
}

VertexData<Vector3> InverseSolver::ARAP_update_positions(VertexData<Vector3> hull_updates){
    size_t n = forwardSolver->inputMesh->nVertices();
    size_t max_iters = 1; // doesnt need much given the small updates
    VertexData<Vector3> old_positions = forwardSolver->inputGeometry->inputVertexPositions,
                        new_positions = forwardSolver->inputGeometry->inputVertexPositions;
    Vector<size_t> hull_indices(forwardSolver->hullMesh->nVertices()), // will be used for row selection later
                   interior_indices(n - forwardSolver->hullMesh->nVertices());

    for (Vertex v: forwardSolver->hullMesh->vertices()){
        Vertex interior_vertex = forwardSolver->inputMesh->vertex(forwardSolver->org_hull_indices[v]);
        new_positions[interior_vertex] = old_positions[interior_vertex] + hull_updates[v];
        hull_indices[v.getIndex()] = interior_vertex.getIndex();
    }
    size_t cnt = 0;
    for (Vertex v: forwardSolver->inputMesh->vertices()){
        if (forwardSolver->on_hull_index[v] == INVALID_IND){
            interior_indices[cnt++] = v.getIndex();
        }
    }

    // get laplacian and prefactor
    forwardSolver->inputGeometry->requireCotanLaplacian();
    SparseMatrix<double> L(forwardSolver->inputGeometry->cotanLaplacian);
    SparseMatrix<double> Ls = L + 1e-7 * identityMatrix<double>(n);
    PositiveDefiniteSolver<double> L_solver(Ls);
    // constraint matrix
    DenseMatrix<double> temp_Ls = Ls.toDense();
    SparseMatrix<double> constraint_Ls = temp_Ls(interior_indices, interior_indices).sparseView();
    PositiveDefiniteSolver<double> constraint_L_solver(constraint_Ls);

    DenseMatrix<double> old_pos_mat = vertex_data_to_matrix(old_positions);
    DenseMatrix<double> new_pos_mat = vertex_data_to_matrix(new_positions);
    
    // initial solution
    auto delta = Ls * old_pos_mat; // Ls * old_pos
    DenseMatrix<double> constraint_vec(n,3); // zeros for non-constrained vertices
    constraint_vec.fill(0.);
    constraint_vec(hull_indices, Eigen::all) = new_pos_mat(hull_indices, Eigen::all);
    constraint_vec = Ls * constraint_vec;
    DenseMatrix<double> init_sol = solve_dense_b(&L_solver, delta - constraint_vec); // laplacian edit process
    init_sol(hull_indices, Eigen::all) = new_pos_mat(hull_indices, Eigen::all);

    // polyscope::registerSurfaceMesh("init sol", init_sol,
    //                                            forwardSolver->inputMesh->getFaceVertexList());
    DenseMatrix<double> tmp_pos_mat = init_sol,
                        tmp_interior_pos_mat = init_sol(interior_indices, Eigen::all);
    for (size_t i = 0; i < max_iters; i++){
        // evaluate R_i's
        VertexData<DenseMatrix<double>> rotations = find_rotations(old_pos_mat, tmp_pos_mat);
        // find positions; p'
        DenseMatrix<double> b(n,3);
        //  - LHS
        b.fill(0.);
        for (Vertex v: forwardSolver->inputMesh->vertices()){
            auto p_v = old_pos_mat.row(v.getIndex());
            DenseMatrix<double> R_v = rotations[v];
            for (Vertex fv: v.adjacentVertices()){
                DenseMatrix<double> R_fv = rotations[fv];
                auto p_fv = old_pos_mat.row(fv.getIndex());
                b(v.getIndex(), Eigen::all) += (-L.coeffRef(v.getIndex(), fv.getIndex()) * 
                                       (R_v + R_fv) * (p_v - p_fv).transpose()).transpose()/2.;
            }
        }
        // solving with constraints
        DenseMatrix<double> constraint_b = (b - constraint_vec)(interior_indices, Eigen::all);
        tmp_interior_pos_mat = solve_dense_b(&constraint_L_solver, constraint_b);
        tmp_pos_mat(interior_indices, Eigen::all) = tmp_interior_pos_mat;
        polyscope::registerSurfaceMesh("z ARAP iter " + std::to_string(i), tmp_pos_mat,
                                               forwardSolver->inputMesh->getFaceVertexList());
    }
    forwardSolver->inputGeometry->unrequireCotanLaplacian();
    VertexData<Vector3> update_vecs(*forwardSolver->inputMesh);
    for (Vertex v: forwardSolver->inputMesh->vertices()){
        update_vecs[v] = vec2vec3(tmp_pos_mat.row(v.getIndex())) - old_positions[v];
    }
    return update_vecs;
}

VertexData<Vector3> InverseSolver::trivial_update_positions(VertexData<Vector3> hull_updates){
    VertexData<Vector3> updates(*forwardSolver->inputMesh, Vector3::zero());
    for (Vertex v: forwardSolver->hullMesh->vertices()){
        Vertex interior_vertex = forwardSolver->inputMesh->vertex(forwardSolver->org_hull_indices[v]);
        updates[interior_vertex] = hull_updates[v];
    }
    return updates;
}

VertexData<Vector3> InverseSolver::greedy_update_positions(VertexData<Vector3> hull_updates){
    VertexData<Vector3> old_positions = forwardSolver->hullGeometry->inputVertexPositions;
    // temp update
    forwardSolver->hullGeometry->inputVertexPositions += hull_updates;
    printf("hererer\n");
    Vector3 O = forwardSolver->get_G(); // not updated yet!
    VertexData<Vector3> greedy_updates(*forwardSolver->inputMesh);
    for (Vertex v: forwardSolver->inputMesh->vertices()){
        if (forwardSolver->on_hull_index[v] == INVALID_IND){ // interior check!
            Face corres_hull_face = interior_v_to_hull_f[v];
            Vector3 corresF_new_normal = forwardSolver->hullGeometry->faceNormal(corres_hull_face),
                    corresF_old_normal = old_normals[corres_hull_face.getIndex()];
            Vector3 old_OP = forwardSolver->inputGeometry->inputVertexPositions[v] - O;
            Vector3 rot_axis_unnormalized = cross(corresF_old_normal, corresF_new_normal);
            double rotation_theta = asin(rot_axis_unnormalized.norm());
            Vector3 new_OP = old_OP.rotateAround(rot_axis_unnormalized.normalize(), rotation_theta);
            
            // possible re-scale; TODO: always rescale even if not hitting
            std::vector<Vector3> polygon_points;
            for (Vertex fv: corres_hull_face.adjacentVertices())
                polygon_points.push_back(forwardSolver->hullGeometry->inputVertexPositions[fv]);
            double new_t = ray_intersect(O, new_OP.normalize(), polygon_points);
            if (new_t != -1.)
                new_OP = new_OP.normalize() * interior_v_to_hull_f_hit_ratio[v] * new_t;
            // o.w. just rotate

            greedy_updates[v] = O + new_OP - old_positions[v];
        }
        else 
            greedy_updates[v] = hull_updates[forwardSolver->on_hull_index[v]];
    }
    // undo hull updates?
    forwardSolver->hullGeometry->inputVertexPositions = old_positions;
    return greedy_updates;
}

VertexData<Vector3> InverseSolver::diffusive_update_positions(VertexData<Vector3> hull_updates){
    //  --- updating with diffusion ---  //
    forwardSolver->inputGeometry->requireEdgeLengths();
    double meanEdgeLength = 0.;
    for (Edge e : forwardSolver->inputMesh->edges()) {
        meanEdgeLength += forwardSolver->inputGeometry->edgeLengths[e];
    }
    meanEdgeLength /= forwardSolver->inputMesh->nEdges();
    double shortTime = 0.1 * meanEdgeLength * meanEdgeLength;

    forwardSolver->inputGeometry->requireVertexLumpedMassMatrix();
    SparseMatrix<double>& M = forwardSolver->inputGeometry->vertexLumpedMassMatrix;

    // Laplacian
    forwardSolver->inputGeometry->requireCotanLaplacian();
    SparseMatrix<double> L = forwardSolver->inputGeometry->cotanLaplacian;

    // Heat operator
    SparseMatrix<double> heatOp = M + shortTime * L;

    PositiveDefiniteSolver<double> heatSolver(heatOp);
    
    // Build RHS
    VertexData<double> rhsx(*forwardSolver->inputMesh, 0.),
                       rhsy(*forwardSolver->inputMesh, 0.),
                       rhsz(*forwardSolver->inputMesh, 0.);
    for (Vertex hull_v: forwardSolver->hullMesh->vertices()){
        Vertex interior_vertex = forwardSolver->inputMesh->vertex(forwardSolver->org_hull_indices[hull_v]); // asserts correct assignment of hull indices
        rhsx[forwardSolver->org_hull_indices[hull_v]] = hull_updates[hull_v].x;
        rhsy[forwardSolver->org_hull_indices[hull_v]] = hull_updates[hull_v].y;
        rhsz[forwardSolver->org_hull_indices[hull_v]] = hull_updates[hull_v].z;
    }
    Vector<double> rhsx_Vec = rhsx.toVector(),
                   rhsy_Vec = rhsy.toVector(),
                   rhsz_Vec = rhsz.toVector();
    Vector<double> diffuse_x_Vec = heatSolver.solve(rhsx_Vec),
                   diffuse_y_Vec = heatSolver.solve(rhsy_Vec),
                   diffuse_z_Vec = heatSolver.solve(rhsz_Vec);

    VertexData<Vector3> diffused_updates(*forwardSolver->inputMesh);
    for (Vertex v: forwardSolver->inputMesh->vertices()){
        diffused_updates[v] = Vector3{diffuse_x_Vec[v.getIndex()], 
                                      diffuse_y_Vec[v.getIndex()], 
                                      diffuse_z_Vec[v.getIndex()]};
    }
    // keep the original updates ?
    for (Vertex hull_v: forwardSolver->hullMesh->vertices()){
        Vertex interior_vertex = forwardSolver->inputMesh->vertex(forwardSolver->org_hull_indices[hull_v]); // asserts correct assignment of hull indices
        diffused_updates[interior_vertex] = hull_updates[hull_v];
    }

    // un-require stuff
    forwardSolver->inputGeometry->unrequireEdgeLengths();
    forwardSolver->inputGeometry->unrequireVertexLumpedMassMatrix();
    forwardSolver->inputGeometry->unrequireCotanLaplacian();

    return diffused_updates;
}