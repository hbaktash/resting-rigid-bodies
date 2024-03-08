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


// optimization stuff

InverseSolver::InverseSolver(BoundaryBuilder* boundaryBuilder){
    this->boundaryBuilder = boundaryBuilder;
    forwardSolver = boundaryBuilder->forward_solver;
    set_fair_distribution();
    save_initial_pos_and_Ls(); // calling required L
    // initialize_constrained_L_solver(forwardSolver->interior_indices); // calling required L
}

void InverseSolver::save_initial_pos_and_Ls(){
    initial_geometry =  new VertexPositionGeometry(*forwardSolver->inputMesh,
                                                   vertex_data_to_matrix(forwardSolver->inputGeometry->inputVertexPositions));
    initial_geometry->requireCotanLaplacian(); // could save this in future
    initial_Ls = initial_geometry->cotanLaplacian;
    initial_geometry->unrequireCotanLaplacian();
    current_Ls = initial_Ls;
}

void InverseSolver::set_fair_distribution() {
    goal_prob = FaceData<double>(*forwardSolver->hullMesh,
                                 1./(double)forwardSolver->hullMesh->nFaces());
}

void InverseSolver::set_fair_distribution_for_sink_faces(size_t goal_stable_count){
    std::vector<double> face_areas;
    size_t count = goal_stable_count;
    goal_prob = FaceData<double>(*forwardSolver->hullMesh, 0.);
    for (Face f: forwardSolver->hullMesh->faces()){
        if (flow_structure[f] == f) { // f is sink
            face_areas.push_back(boundaryBuilder->face_region_area[f]);
        }
    }
    // select the high prob faces for optimizations
    std::sort(face_areas.begin(), face_areas.end());
    double nth_largest_area = face_areas[face_areas.size() - std::min(goal_stable_count, face_areas.size()) ];
    for (Face f: forwardSolver->hullMesh->faces()){
        Vector3 face_normal = forwardSolver->hullGeometry->faceNormal(f);
        if (flow_structure[f] == f && boundaryBuilder->face_region_area[f] >= nth_largest_area){
            goal_prob[f] = 1./(double)count;
            old_stable_normals.push_back(face_normal);
            // printf(" face %d goal area: %f \n", f.getIndex(), goal_prob[f]);
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
        // if (closest_normal_dist < normal_threshold){
            printf("    - updated to face %d\n", closest_face.getIndex());
            new_stable_faces.push_back(closest_face);
            new_stables_count++;
            new_stable_normals.push_back(closest_f_normal);
        // }
        // else {
        //     printf("    - face %d removed due to distance %f\n", closest_face.getIndex(), closest_normal_dist);
        //     polyscope::warning("    - f" + std::to_string(closest_face.getIndex()) + " removed due to dist" + std::to_string(closest_normal_dist)+ "\n");
        // }
    }
    double goal_fair_prob = 1./(double)new_stables_count;
    for (Face f: new_stable_faces){
        goal_prob[f] += goal_fair_prob; // += in case two stables merge into one
        printf(" face %d goal area: %f \n", f.getIndex(), goal_prob[f]);
    }
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
void InverseSolver::find_uni_mass_d_pf_dv(bool frozen_G, bool check_FD){
    //  assuming that we only move convex hull vertices
    //  and that input is convex
    if (!frozen_G){
        find_d_pf_d_Gs();
        find_dG_dvs();
    }
    find_d_pf_dvs();

    uni_mass_d_pf_dv = FaceData<VertexData<Vector3>>(*forwardSolver->hullMesh);
    Vector3 zvec = Vector3::zero();
    Vector3 G = forwardSolver->get_G();
    // printf("here\n");
    for (Face f: forwardSolver->hullMesh->faces()){
        uni_mass_d_pf_dv[f] = VertexData<Vector3>(*forwardSolver->hullMesh, zvec);
        if (!frozen_G && compute_global_G_effect){
            for (Vertex v: forwardSolver->hullMesh->vertices()){
                // center of mass only term
                Vector3 term1 = d_pf_dv[f][v];
                
                Vertex org_vertex = forwardSolver->inputMesh->vertex(forwardSolver->org_hull_indices[v]);
                // geometric term: gauss map geometry
                Vector<double> term2 = dG_dv[org_vertex].transpose() * vec32vec(d_pf_d_G[f]);
                uni_mass_d_pf_dv[f][v] = term1 + Vector3({term2[0], term2[1], term2[2]});
            }
        }
        else {
            for (Vertex v: f.adjacentVertices()){
                if (frozen_G)
                    uni_mass_d_pf_dv[f][v] = d_pf_dv[f][v];
                else {
                    Vertex org_vertex = forwardSolver->inputMesh->vertex(forwardSolver->org_hull_indices[v]);
                    Vector<double> term2 = dG_dv[org_vertex].transpose() * vec32vec(d_pf_d_G[f]);
                    uni_mass_d_pf_dv[f][v] = d_pf_dv[f][v] + Vector3({term2[0], term2[1], term2[2]});
                }
            }
        }
    }
    if (check_FD){
        printf("FD check: \n");
        double step = 1e-6;
        if (!frozen_G)
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
                if (!frozen_G)
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
VertexData<Vector3> InverseSolver::find_uni_mass_total_vertex_grads(size_t goal_stable_count,
                                                                    bool with_flow_structure,                        
                                                                    double stable_normal_update_thres){
    Vector3 zvec = Vector3::zero();
    VertexData<Vector3> total_vertex_grads(*forwardSolver->hullMesh, zvec);
    // if (with_flow_structure && stable_normal_update_thres < 0){
    //     set_fair_distribution_for_sink_faces(); // 
    // }
    // else if (stable_normal_update_thres > 0)
    //     update_fair_distribution(stable_normal_update_thres);
    set_fair_distribution_for_sink_faces(goal_stable_count); // TODO: TESTING

    for (Face f: forwardSolver->hullMesh->faces()){ 
        Face last_sink_face;
        double goal_multiplier = with_flow_structure ? 
                        (goal_prob[flow_structure[f]]* 4.*PI - boundaryBuilder->face_region_area[flow_structure[f]]):
                        (goal_prob[f]* 4.*PI - boundaryBuilder->face_region_area[f]);
        if (compute_global_G_effect) {
            for (Vertex v: forwardSolver->hullMesh->vertices()){
                total_vertex_grads[v] += uni_mass_d_pf_dv[f][v] * 
                                                goal_multiplier;
            }
        }
        else {
            for (Vertex v: f.adjacentVertices()){
                total_vertex_grads[v] += uni_mass_d_pf_dv[f][v] * 
                                                goal_multiplier;
                // printf("      -- tmp grad norm for v %d: %f \n", v.getIndex(), uni_mass_d_pf_dv[f][v].norm());
            }
        }
        // printf("   finding total grad for f %d area: %f goal area: %f \n", f.getIndex(),
        //                                                                 boundaryBuilder->face_region_area[last_sink_face]/(4.*PI),
        //                                                                 goal_prob[last_sink_face]);
    }

    return total_vertex_grads;
}


void InverseSolver::subdivide_for_aggressive_updates(VertexData<Vector3> hull_updates){
    // sub-dividing for aggressive updates
    std::vector<Edge> to_split_edges;
    for (Edge e: forwardSolver->inputMesh->edges()){
        Vertex v1 = e.firstVertex(), 
                v2 = e.secondVertex();
        size_t on_hull_idx1 = forwardSolver->on_hull_index[v1],
               on_hull_idx2 = forwardSolver->on_hull_index[v2];
        if (on_hull_idx1 != INVALID_IND && on_hull_idx2 != INVALID_IND){
            Vertex on_hull_v1 = forwardSolver->hullMesh->vertex(on_hull_idx1),
                   on_hull_v2 = forwardSolver->hullMesh->vertex(on_hull_idx2);
            if (hull_updates[on_hull_v1].norm() != 0. && hull_updates[on_hull_v2].norm() != 0.){
                to_split_edges.push_back(e);
            }
        }
    }
    printf(" splitting %d edges\n", to_split_edges.size());
    printf("current mesh size: %d\n", forwardSolver->inputMesh->nVertices());
    for (Edge e: to_split_edges){
        Vertex v1 = e.firstVertex(), 
               v2 = e.secondVertex();
        Halfedge he = forwardSolver->inputMesh->splitEdgeTriangular(e);
        // only containers that need update should be the index-trackers
        forwardSolver->on_hull_index[he.vertex()] = INVALID_IND; // made sure it lays inside the hull
        // update the current mesh-geometry
        forwardSolver->inputGeometry->inputVertexPositions[he.vertex()] = 
            (forwardSolver->inputGeometry->inputVertexPositions[v1] + 
             forwardSolver->inputGeometry->inputVertexPositions[v2])/2.;
        // update the inital mesh-geometry; to get the initial laplacian
        initial_geometry->inputVertexPositions[he.vertex()] = 
            (initial_geometry->inputVertexPositions[v1] + 
             initial_geometry->inputVertexPositions[v2])/2.;
        Vector3 curr_inward_offset = (forwardSolver->inputGeometry->faceNormal(he.face()) + 
                                      forwardSolver->inputGeometry->faceNormal(he.twin().face())) * 1e-6 * -1.,
                initial_inward_offset = (initial_geometry->faceNormal(he.face()) + 
                                         initial_geometry->faceNormal(he.twin().face())) * 1e-6 * -1.; // assuming outward normals
        forwardSolver->inputGeometry->inputVertexPositions[he.vertex()] += curr_inward_offset;
        initial_geometry->inputVertexPositions[he.vertex()] += initial_inward_offset;
    }
    // resizing laplacian and stuff
    forwardSolver->inputGeometry->refreshQuantities();
    initial_geometry->refreshQuantities();

    // updating interior indices
    forwardSolver->update_hull_index_arrays();
    
    initial_geometry->requireCotanLaplacian();
    initial_Ls = initial_geometry->cotanLaplacian;
    initial_geometry->unrequireCotanLaplacian();
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
    for (Vertex v: forwardSolver->inputMesh->vertices()){
        size_t n_v = v.degree();
        Vector<size_t> neigh_inds(n_v);
        Vector<double> weights(n_v);
        size_t cnt = 0;
        for (Edge e: v.adjacentEdges()){
            Vertex neigh_v = e.otherVertex(v);
            neigh_inds(cnt) = neigh_v.getIndex();
            weights(cnt) = current_Ls.coeff(e.firstVertex().getIndex(), e.secondVertex().getIndex());
            cnt++;
        }
        assert(n_v == cnt);
        DenseMatrix<double> weights_diag = weights.asDiagonal();
        DenseMatrix<double> P_v = old_pos(v.getIndex(), Eigen::all).replicate(n_v, 1) - 
                                   old_pos(neigh_inds, Eigen::all),
                            P_v_new = new_pos(v.getIndex(), Eigen::all).replicate(n_v, 1) - 
                                       new_pos(neigh_inds, Eigen::all);
        DenseMatrix<double> S_v = P_v.transpose() * weights_diag * P_v_new;

        auto bdcSVD = S_v.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
        DenseMatrix<double> U = bdcSVD.matrixU(),
                            V = bdcSVD.matrixV();
        DenseMatrix<double> R_v = V * U.transpose();
        rotations[v] = R_v;
        // std::cout<< "rot at v " << v.getIndex() << ":\n" << R_v << "\n";
    }
    return rotations;
}

void InverseSolver::update_Ls(){
    size_t n = forwardSolver->inputMesh->nVertices();
    printf("old geometry flag: %d\n", use_old_geometry);
    if (use_old_geometry)
        current_Ls = initial_Ls;
    else {
        printf("recomputing Ls\n");
        // forwardSolver->inputGeometry->unrequireCotanLaplacian();
        forwardSolver->inputGeometry->requireCotanLaplacian();
        current_Ls = forwardSolver->inputGeometry->cotanLaplacian + 1e-7 * identityMatrix<double>(n);
        forwardSolver->inputGeometry->unrequireCotanLaplacian();
    }
}


DenseMatrix<double> InverseSolver::solve_constrained_Laplace(Vector<size_t> interior_indices, 
                                                             DenseMatrix<double> old_pos, DenseMatrix<double> new_pos,
                                                             bool update_solver_decomp){
    size_t n = forwardSolver->inputMesh->nVertices();
    // get laplacian and prefactor
    geometrycentral::SparseMatrix<double> Ls = current_Ls;
    
    // constrained matrix
    if (update_solver_decomp){
        Vector<bool> interior_indicator(n);
        interior_indicator.fill(false);
        printf("making interior indicator\n");
        for (size_t int_idx: interior_indices)
            interior_indicator[int_idx] = true;
        printf("interior indicator done! Ls size %d %d, indic size %d\n", Ls.rows(), Ls.cols(), interior_indicator.size());
        BlockDecompositionResult<double> decomp = blockDecomposeSquare(Ls, interior_indicator, false);
        printf("decomposed!\n");
        constrained_L_solver = new PositiveDefiniteSolver<double>(decomp.AA);
        printf("solver created!\n");
    }

    // initial solution
    DenseMatrix<double> padded_handles(n,3); // zeros for un-constrained vertices
    padded_handles.fill(0.);
    padded_handles(forwardSolver->hull_indices, Eigen::all) = new_pos(forwardSolver->hull_indices, Eigen::all);
    DenseMatrix<double> sol(n,3);
    // update constraint solver since hull indices changed
    sol(interior_indices, Eigen::all) = solve_dense_b(constrained_L_solver, (Ls*(old_pos - padded_handles))(interior_indices, Eigen::all)); // laplacian edit process
    sol(forwardSolver->hull_indices, Eigen::all) = new_pos(forwardSolver->hull_indices, Eigen::all);
    return sol;
}


VertexData<Vector3> InverseSolver::laplace_update_positions(VertexData<Vector3> new_positions){

    VertexData<Vector3> old_positions = forwardSolver->inputGeometry->inputVertexPositions;
    printf("solving with dirichlet condiitons\n");
    if (use_old_geometry)
        old_positions = initial_geometry->inputVertexPositions;
    // could also use initial positions, but unclear how to index things since convex vectices are not the same
    DenseMatrix<double> old_pos_mat = vertex_data_to_matrix(old_positions); 
    DenseMatrix<double> new_pos_mat = vertex_data_to_matrix(new_positions);
    
    update_Ls(); // whether to use initial Laplacian or the current one
    DenseMatrix<double> sol = solve_constrained_Laplace(forwardSolver->interior_indices, old_pos_mat, new_pos_mat);
    
    VertexData<Vector3> update_vecs(*forwardSolver->inputMesh);
    for (Vertex v: forwardSolver->inputMesh->vertices()){
        update_vecs[v] = vec2vec3(sol.row(v.getIndex())) - old_positions[v];
    }
    return update_vecs;
}


VertexData<Vector3> InverseSolver::ARAP_update_positions(VertexData<Vector3> new_positions){
    size_t n = forwardSolver->inputMesh->nVertices();
    size_t max_iters = arap_max_iter; // shouldn't need much given the small updates?

    VertexData<Vector3> old_positions = forwardSolver->inputGeometry->inputVertexPositions;
    if (use_old_geometry){
        old_positions = initial_geometry->inputVertexPositions;
        printf("using old positions\n");
    }
    DenseMatrix<double> old_pos_mat = vertex_data_to_matrix(old_positions);
    DenseMatrix<double> new_pos_mat = vertex_data_to_matrix(new_positions);
    
    // constrained matrix
    Vector<size_t> interior_indices = forwardSolver->interior_indices;
    update_Ls(); // whether to use initial Laplacian or the current one
    DenseMatrix<double> init_sol = solve_constrained_Laplace(interior_indices, old_pos_mat, new_pos_mat);;
    
    polyscope::registerSurfaceMesh("z ARAP Laplace init", init_sol,
                                           forwardSolver->inputMesh->getFaceVertexList());
    DenseMatrix<double> tmp_pos_mat = init_sol;
    if (libigl_ARAP){
        tmp_pos_mat = get_ARAP_positions(old_pos_mat, new_pos_mat, init_sol, *forwardSolver->inputMesh, forwardSolver->hull_indices.cast<int>());
    }
    else {
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
                // if (forwardSolver->on_hull_index[v] != INVALID_IND)
                //     R_v = identityMatrix<double>(3).toDense();
                for (Vertex v_neigh: v.adjacentVertices()){
                    DenseMatrix<double> R_v_neigh = rotations[v_neigh];
                    // if (forwardSolver->on_hull_index[v_neigh] != INVALID_IND)
                    //     R_v_neigh = identityMatrix<double>(3).toDense();
                    auto p_v_neigh = old_pos_mat.row(v_neigh.getIndex());
                    b(v.getIndex(), Eigen::all) += current_Ls.coeff(v.getIndex(), v_neigh.getIndex()) * 
                                        ((R_v + R_v_neigh).transpose() * (p_v - p_v_neigh).transpose()).transpose();
                }
            }
            b = b/2.;
            // solving with constraints
            DenseMatrix<double> padded_handles(n,3); // zeros for un-constrained vertices
            padded_handles.fill(0.);
            padded_handles(forwardSolver->hull_indices, Eigen::all) = new_pos_mat(forwardSolver->hull_indices, Eigen::all);
            tmp_pos_mat(interior_indices, Eigen::all) = solve_dense_b(constrained_L_solver,  
                                                (b - current_Ls*padded_handles)(interior_indices, Eigen::all));
            polyscope::registerSurfaceMesh("z ARAP iter " + std::to_string(i), tmp_pos_mat,
                                                forwardSolver->inputMesh->getFaceVertexList())->setEnabled(false);
        }
    }
    VertexData<Vector3> update_vecs(*forwardSolver->inputMesh);
    for (Vertex v: forwardSolver->inputMesh->vertices()){
        update_vecs[v] = vec2vec3(tmp_pos_mat.row(v.getIndex())) - forwardSolver->inputGeometry->inputVertexPositions[v];
    }
    return update_vecs;
}

VertexData<Vector3> InverseSolver::trivial_update_positions(VertexData<Vector3> hull_updates){
    VertexData<Vector3> updates(*forwardSolver->inputMesh, Vector3::zero());
    for (Vertex v: forwardSolver->hullMesh->vertices()){
        Vertex original_vertex = forwardSolver->inputMesh->vertex(forwardSolver->org_hull_indices[v]);
        updates[original_vertex] = hull_updates[v];
    }
    return updates;
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
    geometrycentral::SparseMatrix<double>& M = forwardSolver->inputGeometry->vertexLumpedMassMatrix;

    // Laplacian
    forwardSolver->inputGeometry->requireCotanLaplacian();
    geometrycentral::SparseMatrix<double> L = forwardSolver->inputGeometry->cotanLaplacian;

    // Heat operator
    geometrycentral::SparseMatrix<double> heatOp = M + shortTime * L;

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


VertexData<Vector3> InverseSolver::sobolev_diffuse_gradients(VertexData<Vector3> grads, double sobolev_lambda, size_t sobolev_p){
    size_t n = forwardSolver->hullMesh->nVertices(), 
           e = forwardSolver->hullMesh->nEdges();
    SparseMatrix<double> graph_L(n, n);
    std::vector<Eigen::Triplet<double>> gL_tripletList;
    gL_tripletList.reserve(4*e);
    for (Edge e: forwardSolver->hullMesh->edges()){
        Vertex v1 = e.halfedge().vertex(),
               v2 = e.halfedge().twin().vertex();
        gL_tripletList.push_back(Eigen::Triplet<double>(v1.getIndex(), v2.getIndex(), -1.));
        gL_tripletList.push_back(Eigen::Triplet<double>(v2.getIndex(), v1.getIndex(), -1.));
        gL_tripletList.push_back(Eigen::Triplet<double>(v1.getIndex(), v1.getIndex(), 1.));
        gL_tripletList.push_back(Eigen::Triplet<double>(v2.getIndex(), v2.getIndex(), 1.));
    }
    graph_L.setFromTriplets(gL_tripletList.begin(), gL_tripletList.end());

    // Sobolev operator
    SparseMatrix<double> sobolevOp = graph_L + sobolev_lambda * identityMatrix<double>(n);
    PositiveDefiniteSolver<double> sobolevSolver(sobolevOp);
    DenseMatrix<double> sobolev_grads;
    for (size_t i = 0; i < sobolev_p; i++){
        sobolev_grads = solve_dense_b(&sobolevSolver, vertex_data_to_matrix(grads));
    }
    VertexData<Vector3> sobolev_grads_vd(*forwardSolver->hullMesh);
    for (Vertex v: forwardSolver->hullMesh->vertices()){
        sobolev_grads_vd[v] = vec2vec3(sobolev_grads.row(v.getIndex()));
    }
    return sobolev_grads_vd;
}