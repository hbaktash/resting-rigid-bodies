// Forward declarations


template <typename Scalar>
Scalar BoundaryBuilder::dice_energy(Eigen::MatrixX3<Scalar> hull_positions, Eigen::Vector3<Scalar> G,
                                    Forward3DSolver &tmp_solver, double bary_reg, double co_planar_reg, double cluster_distance_reg,
                                    std::string policy, FaceData<double> goal_probs, 
                                    size_t side_count, bool verbose
                                    ){
    // Eigen::MatrixX3d hull_positions_d = static_cast<Eigen::MatrixX3d>(hull_positions);//.template cast<double>();
    // Eigen::Vector3d G_d = G.template cast<double>();
    // Forward3DSolver tmp_solver(hull_positions_d, G_d, true); // assuming input is convex; will be asserted internally in the constructor
    // // precomputes
    tmp_solver.initialize_pre_computes();
    FaceData<Face> face_last_face = tmp_solver.face_last_face;
    VertexData<bool> vertex_is_stabilizable = tmp_solver.vertex_is_stabilizable; 
    EdgeData<Vertex> edge_next_vertex = tmp_solver.edge_next_vertex;
    // tmp_solver.print_precomputes();
    std::vector<Edge> terminal_edges;
    // printf("  finding terminal edges \n");
    size_t cnt = 0;
    for (Edge e: tmp_solver.hullMesh->edges()) {
        if (tmp_solver.edge_next_vertex[e].getIndex() == INVALID_IND){ // singular edge
            Face f1 = e.halfedge().face(),
                 f2 = e.halfedge().twin().face();
            if (tmp_solver.face_last_face[f1] != tmp_solver.face_last_face[f2]){ // saddle edge
                terminal_edges.push_back(e);
                // printf("cnt %d: edge %d is terminal with faces %d, %d\n", cnt, e.getIndex(), f1.getIndex(), f2.getIndex());
                // printf("    --- face last faces: %d, %d\n", tmp_solver.face_last_face[f1].getIndex(), tmp_solver.face_last_face[f2].getIndex());
                
                cnt++;
            }
        }
    }

    // pre compute face nomrals; GC normals not templated
    FaceData<Eigen::Vector3<Scalar>> face_normals(*tmp_solver.hullMesh);
    for (Face f: tmp_solver.hullMesh->faces()){
        Halfedge he = f.halfedge();
        Vertex v1 = he.vertex(),
               v2 = he.next().vertex(),
               v3 = he.next().next().vertex();
        face_normals[f] = (hull_positions.row(v2.getIndex()) - hull_positions.row(v1.getIndex()))
                         .cross(hull_positions.row(v3.getIndex()) - hull_positions.row(v1.getIndex())).normalized();//.normalized();
    }

    // morse complex areas; only non-zero for stable faces
    FaceData<Scalar> face_region_area(*tmp_solver.hullMesh, 0.);
    
    // for barycentric energy
    FaceData<std::set<Vertex>> face_region_boundary_vertices(*tmp_solver.hullMesh);
    FaceData<Scalar> distance_to_barycenters(*tmp_solver.hullMesh, 0.);
    // printf("  stan: back-flowing terminal edges %lu \n", terminal_edges.size());
    int i = 1;
    for (Edge e: terminal_edges){
        // printf("\n - starting at terminal edge: %d/%d \n", i++, terminal_edges.size());
        // std::cout << "  -e" << e.getIndex() << " : " << e.firstVertex().getIndex() << " , " << e.secondVertex().getIndex() << " adj faces: "<< e.halfedge().face().getIndex() << " " << e.halfedge().twin().face().getIndex() << std::endl;
        // bnd_normal->normal_ad = point_to_segment_normal_ad(e);
        Eigen::Vector3<Scalar> edge_bnd_normal = point_to_segment_normal<Scalar>(G, hull_positions.row(e.firstVertex().getIndex()), hull_positions.row(e.secondVertex().getIndex()));
        
        //
        // std::vector<std::array<size_t, 2>> vertex_pairs;
        // vertex_pairs.push_back({e.firstVertex().getIndex(), e.secondVertex().getIndex()});
        // polyscope::registerCurveNetwork("current edge", hull_positions, vertex_pairs);
        
        Edge current_e = e;
        Face f1 = current_e.halfedge().face(),
             f2 = current_e.halfedge().twin().face();
        
        // //
        // draw_arc_on_sphere_static_temp(tmp_solver.hullGeometry->faceNormal(f1), 
        //                           tmp_solver.hullGeometry->faceNormal(f2), 
        //                           Vector3({0,2,0}), 1., 12, 1, 0.1, glm::vec3(1., 0., 1.), 0.5);
        
        Face ff1 = face_last_face[f1];
        Face ff2 = face_last_face[f2];
        // for visuals
        Eigen::Vector3<Scalar> ff1_normal  = face_normals[ff1], // static
                               ff2_normal  = face_normals[ff2]; // static
        Eigen::Vector3<Scalar> adj_f1_normal  = face_normals[f1], // static
                               adj_f2_normal  = face_normals[f2]; // static
        // std::cout << " ff1: " << ff1.getIndex() << " ff2: " << ff2.getIndex() << std::endl;
        // std::cout << " with normals: " << ff1_normal.transpose() << "  --  " << ff2_normal.transpose() << std::endl;
        for (Vertex v: {e.firstVertex(), e.secondVertex()}){ // two ways to go from a saddle edge
            
            Eigen::Vector3<Scalar> tmp_normal  = edge_bnd_normal;   // changes in the while loop
            Vertex current_v = v;                            // changes in the while loop
            Eigen::Vector3<Scalar> v_normal    = (hull_positions.row(current_v.getIndex()) - G.transpose()).normalized(); // changes in the while loop
            
            // polyscope::registerPointCloud("current V", std::vector<Vector3>{vec_to_GC_vec3(tmp_normal) +Vector3({0,2,0})});
            // polyscope::show();

            // ** DONT CONVERT these to template**
            double ff1_area_sign = ff1_normal.dot(tmp_normal.cross(v_normal)) >= 0 ? 1. : -1.; // ff1 on rhs of bndN->vN; static
            double ff2_area_sign = ff2_normal.dot(tmp_normal.cross(v_normal)) >= 0 ? 1. : -1.; // ff2 on rhs of bndN->vN; static
            double adj_f1_area_sign = adj_f1_normal.dot(tmp_normal.cross(v_normal)) >= 0 ? 1. : -1.; // ff1 on rhs of bndN->vN; static
            double adj_f2_area_sign = adj_f2_normal.dot(tmp_normal.cross(v_normal)) >= 0 ? 1. : -1.; // ff2 on rhs of bndN->vN; static
            
            ff1_area_sign *= adj_f1_area_sign == ff1_area_sign ? 1. : -1.;
            ff2_area_sign *= adj_f2_area_sign == ff2_area_sign ? 1. : -1.;

            Eigen::Vector3<Scalar> next_normal;                     // changes in the while loop
            int cnt = 0;
            while (true) {
                cnt++;
                // std::cout << "   - at vertex " << current_v.getIndex() << std::endl;
                // if (cnt > 50){
                //     std::cout << "  - too many iterations\n";
                //     break;
                // }
                if (vertex_is_stabilizable[current_v]){
                    next_normal = v_normal;
                    Scalar patch_area_ff1 = triangle_patch_signed_area_on_sphere<Scalar>(ff1_normal, tmp_normal, next_normal);
                    face_region_area[ff1] += ff1_area_sign * patch_area_ff1;
                    Scalar patch_area_ff2 = triangle_patch_signed_area_on_sphere<Scalar>(ff2_normal, tmp_normal, next_normal);
                    face_region_area[ff2] += ff2_area_sign * patch_area_ff2;
                    // draw_arc_on_sphere_static(vec2vec3(tmp_normal), vec2vec3(next_normal), Vector3({0,2,0}), 1., 12, 0, 0.1, glm::vec3(1., 0., 1.), 0.5);
                    // polyscope::show();         
                    // for barycentric energy
                    face_region_boundary_vertices[ff1].insert(current_v);
                    face_region_boundary_vertices[ff2].insert(current_v);
                    break; // Exit while; got to a maximum
                }
                else{
                    for (Edge neigh_e: current_v.adjacentEdges()){
                        // printf(" -e %d-", neigh_e.getIndex());
                        if (neigh_e != current_e && edge_next_vertex[neigh_e] == current_v){ // neigh_e is a source for this vertex
                            Face f1 = neigh_e.halfedge().face(),
                                 f2 = neigh_e.halfedge().twin().face();
                            bool sign_change = false;
                            next_normal = intersect_arcs(v_normal, tmp_normal, 
                                                         face_normals[f1], 
                                                         face_normals[f2]);
                            if(next_normal.norm() != 0.){ // found the next source normal
                                Scalar patch_area_ff1 = triangle_patch_signed_area_on_sphere<Scalar>(ff1_normal, tmp_normal, next_normal);
                                face_region_area[ff1] += ff1_area_sign * patch_area_ff1;
                                Scalar patch_area_ff2 = triangle_patch_signed_area_on_sphere<Scalar>(ff2_normal, tmp_normal, next_normal);
                                face_region_area[ff2] += ff2_area_sign * patch_area_ff2;

                                // //
                                // draw_arc_on_sphere_static(vec_to_GC_vec3(tmp_normal), vec_to_GC_vec3(next_normal), Vector3({0,2,0}), 1., 12, 0, 0.1, glm::vec3(1., 0., 1.), 0.5);
                                // polyscope::show();
                                // updates
                                tmp_normal = next_normal;
                                current_v = neigh_e.otherVertex(current_v);
                                current_e = neigh_e;
                                v_normal  = (hull_positions.row(current_v.getIndex()) - G.transpose()).normalized();
                                break; // don't check other edges; go for next cell
                            }
                        }
                    }
                }
            }
        }
    }
    

    //// regularizers
    // compute distance to bary centers
    Scalar bary_energy = 0.;
    Scalar co_planar_energy = 0.;
    // split policy
    std::string policy_shape = policy.substr(policy.find(" ") + 1), // second word
                policy_general = policy.substr(0, policy.find(" ")); // first word
    if (policy_general == "fair" || policy_general == "manual"){
        for (Face f: tmp_solver.hullMesh->faces()){
            if (face_region_area[f] > 0){
                Eigen::Vector3<Scalar> barycenter = Eigen::Vector3<Scalar>::Zero();
                for (Vertex v: face_region_boundary_vertices[f]){
                    barycenter += (hull_positions.row(v.getIndex())- G.transpose()).normalized(); // vertex-G normals
                }
                barycenter /= face_region_boundary_vertices[f].size();
                barycenter = barycenter.normalized();
                distance_to_barycenters[f] = (barycenter - face_normals[f]).norm();
                bary_energy += bary_reg * distance_to_barycenters[f] * distance_to_barycenters[f];
            }
        }
        // penalize non-coplanar faces in the rolling DAG
        for (Face f: tmp_solver.hullMesh->faces()){
            if (face_last_face[f] != f){
                Eigen::Vector3<Scalar> f_normal = face_normals[f];
                Eigen::Vector3<Scalar> last_f_normal = face_normals[face_last_face[f]];
                Scalar normal_diff = (f_normal- last_f_normal).squaredNorm();
                co_planar_energy += co_planar_reg * normal_diff;
            }
        }
    }

    Scalar energy = 0.;
    energy += bary_energy + co_planar_energy;
    if (policy_general == "fair"){
        // sort for dice energy
        Scalar total_area = 0.;
        std::vector<std::pair<Face, Scalar>> probs;
        for (Face f: tmp_solver.hullMesh->faces()){
            if (face_region_area[f] > 0){
                // face_region_area[f] /= (4.*PI);
                probs.push_back({f, face_region_area[f]});
                total_area += face_region_area[f];
            }
        }
        // std::cout << "sum over probabilities: " << total_area/(4.*PI) << std::endl;
        std::sort(probs.begin(), probs.end(), [] (auto a, auto b) { return a.second > b.second; });
        // printf(" static sorted probs: \n");
        size_t tmp_side_cnt = 0;
        double goal_prob = 1. / (double)side_count;
        for (auto pair: probs){ 
            if(verbose){
                std::cout << "  -f" << pair.first.getIndex() 
                          << "(" 
                          << pair.first.halfedge().vertex().getIndex() << ", "
                          << pair.first.halfedge().tipVertex().getIndex()<< ", "
                          << pair.first.halfedge().next().tipVertex().getIndex()
                          << ")" << " : " << pair.second/(4. * PI) << std::endl; 
            }
            if (tmp_side_cnt < side_count){ // goal 1/n
                Scalar diff = face_region_area[pair.first] - goal_prob * 4. * PI;
                energy += diff * diff;
            }
            else { // goal zero
                Scalar diff = face_region_area[pair.first];
                energy += diff * diff;
            }
            tmp_side_cnt++;
        }
        if (tmp_side_cnt < side_count){
            printf("  - not enough faces\n");
            energy += (side_count - tmp_side_cnt) * goal_prob * goal_prob * 4. * PI * 4. * PI;
        }

        // return probs[0].second; // DEBUG: max prob
    }
    else if (policy_general == "manual"){
        FaceData<double> my_probs = manual_stable_only_face_prob_assignment(&tmp_solver, policy_shape);
        for (Face f: tmp_solver.hullMesh->faces()){
            // // norm 2
            Scalar diff = face_region_area[f] - my_probs[f] * 4. * PI;
            energy += diff * diff;
            // // KL divergence
            //     double goal_area = my_probs[f] * 4. * PI;
            //     energy += goal_area * log(goal_area / face_region_area[f]);
            // }
        }
        if (verbose){
            std::cout << "  - bary energy: " << bary_energy << "\t co-planar energy: " << co_planar_energy << std::endl;
            std::cout << "  - pure dice energy: " << energy - bary_energy - co_planar_energy << std::endl;
        }
    }
    else if (policy_general == "manualCluster"){
        std::vector<std::pair<std::vector<Face>, double>> clustered_probs = manual_clustered_face_prob_assignment(&tmp_solver, policy_shape);
        Scalar total_coplanar_energy = 0.,
               pure_dice_energy = 0.,
               total_bary_energy = 0.,
               total_cluster_distance_energy = 0.;
        std::vector<Eigen::Vector3<Scalar>> cluster_barycenters;

        size_t non_zero_clusters = 0;
        for (auto cluster: clustered_probs){

            if (cluster.first.size() > 0)
                non_zero_clusters++;
            else
                continue;
            // Pure dice energy
            Scalar cluster_area_sum = 0.;
            for (Face f: cluster.first){
                if (face_region_area[f] > 0){
                    // sum of stable face areas in the cluster
                    cluster_area_sum += face_region_area[f];
                }
            }
            Scalar diff = cluster_area_sum - cluster.second * 4. * PI;
            pure_dice_energy += diff * diff;

            // Co-planar energy
            // TODO either stables only or all faces in the cluster
            Scalar cluster_coplanar_energy = 0.;
            std::vector<Face> cluster_faces = cluster.first;
            for (Face f1: cluster_faces){
                Face ff1 = face_last_face[f1];
                if (ff1 == f1) // Could remove
                    if (std::find(cluster_faces.begin(), cluster_faces.end(), ff1) != cluster_faces.end()){ // final face is in the same cluster
                        for (Face f2: cluster_faces){
                            Face ff2 = face_last_face[f2];
                            if (ff2 == f2) // Could remove
                                if (std::find(cluster_faces.begin(), cluster_faces.end(), ff2) != cluster_faces.end()){ // final face is in the same cluster
                                    Scalar normal_diff = (face_normals[f1] - face_normals[f2]).squaredNorm();
                                    cluster_coplanar_energy += normal_diff;
                                }
                        }
                    }
            }
            total_coplanar_energy += cluster_coplanar_energy;
            
            // Bary energy
            std::set<Vertex> cluster_boundary_vertices;
            Eigen::Vector3<Scalar> stable_faces_cluster_barycenter = Eigen::Vector3<Scalar>::Zero();
            double stable_face_count = 0;
            if (cluster.first.size() == 0){ // avoiding nan
                continue;
            }
            for (Face cluster_f: cluster.first){
                if (face_last_face[cluster_f] == cluster_f){ // stable cluster faces
                    cluster_boundary_vertices.insert(face_region_boundary_vertices[cluster_f].begin(), face_region_boundary_vertices[cluster_f].end());
                    stable_faces_cluster_barycenter += face_normals[cluster_f];
                    stable_face_count += 1.;
                }
            }
            if (stable_face_count == 0){ // avoiding nan
                continue;
            }
            stable_faces_cluster_barycenter = stable_faces_cluster_barycenter.normalized();
            Eigen::Vector3<Scalar> cluster_barycenter = Eigen::Vector3<Scalar>::Zero();
            for (Vertex v: cluster_boundary_vertices)
                cluster_barycenter += (hull_positions.row(v.getIndex()) - G.transpose()).normalized(); // vertex-G normals
            cluster_barycenter /= (double)cluster_boundary_vertices.size();
            cluster_barycenter = cluster_barycenter.normalized();
            Scalar stable_face_distances_sum_to_barycenter = 0.;
            // -- Penalize very stable face distance to bary center; could de-motivate splitting; // sum (fi - bary)^2
            // for (Face cluster_f: cluster.first)
            //     if (face_region_area[cluster_f] > 0) // stable cluster faces
            //         stable_face_distances_sum_to_barycenter += (cluster_barycenter - face_normals[cluster_f]).norm();
            // -- Penalize only stable barycenters distance to cluster barycenter; might leave splits alone // [mean (fi) - bary]^2
            stable_face_distances_sum_to_barycenter = (cluster_barycenter - stable_faces_cluster_barycenter).norm();
            total_bary_energy += stable_face_distances_sum_to_barycenter * stable_face_distances_sum_to_barycenter;
            
            // cluster distance term
            cluster_barycenters.push_back(stable_faces_cluster_barycenter);
        }
        // cluster distance term
        for (int i = 0; i < cluster_barycenters.size(); i++){
            for (int j = i+1; j < cluster_barycenters.size(); j++){
                Scalar diff = (cluster_barycenters[i] - cluster_barycenters[j]).norm();
                total_cluster_distance_energy += 1./(diff * diff);
            }
        }
        // total energy
        energy += pure_dice_energy + 
                  co_planar_reg * total_coplanar_energy + bary_reg * total_bary_energy + cluster_distance_reg * total_cluster_distance_energy;
        if (verbose){
            // for (int i = 0; i < cluster_barycenters.size(); i++){
            //     std::cout << "  - cluster " << i << " barycenter: " << cluster_barycenters[i].transpose() << std::endl;
            //     if (i < cluster_barycenters.size())
            //         std::cout << "cluster goal prob: " << clustered_probs[i].second << std::endl;
            // }
            std::cout << "  - barycenter dist energy : " << bary_reg << " * " <<  total_bary_energy << std::endl;
            std::cout << "  - coplanar energy        : " << co_planar_reg << " * " << total_coplanar_energy << std::endl;
            std::cout << "  - cluster distance energy: " << cluster_distance_reg << " * " << total_cluster_distance_energy << std::endl;
            std::cout << "  - pure dice energy: " << pure_dice_energy << "\t Nz cluster count: " << non_zero_clusters << std::endl;
            std::cout << "\t\t\t** total energy: " << energy << std::endl;
        }
    }
    else {
        throw std::invalid_argument("policy not recognized");
    }
    return energy;
}

template <typename Scalar>
Eigen::Vector3<Scalar> BoundaryBuilder::point_to_segment_normal(Eigen::Vector3<Scalar> P, Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B){
    Eigen::Vector3<Scalar> PB = B - P,
                    AB = B - A;
    return (PB * AB.dot(AB) - AB * AB.dot(PB)).normalized();
}

template <typename Scalar>
Eigen::Vector3<Scalar> BoundaryBuilder::intersect_arcs(Eigen::Vector3<Scalar> v_normal, Eigen::Vector3<Scalar> R2, Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B){
    Eigen::Vector3<Scalar> p = (((v_normal.cross(R2)).cross(A.cross(B)))).normalized(); // var_G.transpose()
    if (p.cross(B).dot(A.cross(B)) >= 0. && p.cross(A).dot(B.cross(A)) >= 0.)
        return p;
    else if ((-p).cross(B).dot(A.cross(B)) >= 0. && (-p).cross(A).dot(B.cross(A)) >= 0.)
        return -p;
    else
        return Eigen::Vector3<Scalar>::Zero();
    // if (p.dot(A) >= A.dot(B) && p.dot(B) >= A.dot(B))
    //     return p;
    // else if (-p.dot(A) >= A.dot(B) && -p.dot(B) >= A.dot(B))
    //     return -p;
    // else 
    //     return Eigen::Vector3<Scalar>::Zero();
}

template <typename Scalar>
Scalar BoundaryBuilder::triangle_patch_signed_area_on_sphere(Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B, Eigen::Vector3<Scalar> C){
    Scalar Adotbc = A.dot(B.cross(C));
    Scalar denom = A.norm()*B.norm()*C.norm() + A.dot(B)*C.norm() + B.dot(C)*A.norm() + C.dot(A)*B.norm();
    return 2. * atan2(Adotbc, denom);
}
