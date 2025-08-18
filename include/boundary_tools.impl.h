
// #ifdef RESTING_RIGID_BODIES_ENABLE_INVERSE_DESIGN

// template <typename Scalar>
// Scalar BoundaryBuilder::dice_energy(Eigen::MatrixX3<Scalar> hull_positions, Eigen::Vector3<Scalar> G,
//                                     Forward3DSolver &tmp_solver, double bary_reg, double co_planar_reg, double cluster_distance_reg, double unstable_attaction_thresh,
//                                     std::string policy_general, std::vector<std::pair<Vector3, double>> normal_prob_assignment, 
//                                     size_t side_count, bool verbose
//                                     ){
//     // Eigen::MatrixX3d hull_positions_d = static_cast<Eigen::MatrixX3d>(hull_positions);//.template cast<double>();
//     // Eigen::Vector3d G_d = G.template cast<double>();
//     // Forward3DSolver tmp_solver(hull_positions_d, G_d, true); // assuming input is convex; will be asserted internally in the constructor
//     // // precomputes
//     tmp_solver.initialize_pre_computes();
//     FaceData<Face> face_last_face = tmp_solver.face_last_face;
//     VertexData<bool> vertex_is_stabilizable = tmp_solver.vertex_is_stabilizable; 
//     EdgeData<Vertex> edge_next_vertex = tmp_solver.edge_next_vertex;
//     // tmp_solver.print_precomputes();
//     std::vector<Edge> terminal_edges;
//     // printf("  finding terminal edges \n");
//     size_t cnt = 0;
//     for (Edge e: tmp_solver.hullMesh->edges()) {
//         if (tmp_solver.edge_next_vertex[e].getIndex() == INVALID_IND){ // singular edge
//             Face f1 = e.halfedge().face(),
//                  f2 = e.halfedge().twin().face();
//             if (tmp_solver.face_last_face[f1] != tmp_solver.face_last_face[f2]){ // saddle edge
//                 terminal_edges.push_back(e);
//                 // printf("cnt %d: edge %d is terminal with faces %d, %d\n", cnt, e.getIndex(), f1.getIndex(), f2.getIndex());
//                 // printf("    --- face last faces: %d, %d\n", tmp_solver.face_last_face[f1].getIndex(), tmp_solver.face_last_face[f2].getIndex());
                
//                 cnt++;
//             }
//         }
//     }

//     // pre compute face nomrals; GC normals not templated
//     FaceData<Eigen::Vector3<Scalar>> face_normals(*tmp_solver.hullMesh);
//     for (Face f: tmp_solver.hullMesh->faces()){
//         Halfedge he = f.halfedge();
//         Vertex v1 = he.vertex(),
//                v2 = he.next().vertex(),
//                v3 = he.next().next().vertex();
//         face_normals[f] = (hull_positions.row(v2.getIndex()) - hull_positions.row(v1.getIndex()))
//                          .cross(hull_positions.row(v3.getIndex()) - hull_positions.row(v1.getIndex())).normalized();//.normalized();
//     }

//     // morse complex areas; only non-zero for stable faces
//     FaceData<Scalar> face_region_area(*tmp_solver.hullMesh, 0.);
    
//     // for barycentric energy
//     FaceData<std::set<Vertex>> face_region_boundary_vertices(*tmp_solver.hullMesh);
//     FaceData<Scalar> distance_to_barycenters(*tmp_solver.hullMesh, 0.);
//     // printf("  stan: back-flowing terminal edges %lu \n", terminal_edges.size());
//     int i = 1;
//     for (Edge e: terminal_edges){
//         // printf("\n - starting at terminal edge: %d/%d \n", i++, terminal_edges.size());
//         // std::cout << "  -e" << e.getIndex() << " : " << e.firstVertex().getIndex() << " , " << e.secondVertex().getIndex() << " adj faces: "<< e.halfedge().face().getIndex() << " " << e.halfedge().twin().face().getIndex() << std::endl;
//         // bnd_normal->normal_ad = point_to_segment_normal_ad(e);
//         Eigen::Vector3<Scalar> edge_bnd_normal = point_to_segment_normal<Scalar>(G, hull_positions.row(e.firstVertex().getIndex()), hull_positions.row(e.secondVertex().getIndex()));
        
//         //
//         // std::vector<std::array<size_t, 2>> vertex_pairs;
//         // vertex_pairs.push_back({e.firstVertex().getIndex(), e.secondVertex().getIndex()});
//         // polyscope::registerCurveNetwork("current edge", hull_positions, vertex_pairs);
        
//         Edge current_e = e;
//         Face f1 = current_e.halfedge().face(),
//              f2 = current_e.halfedge().twin().face();
        
//         // //
//         // draw_arc_on_sphere_static_temp(tmp_solver.hullGeometry->faceNormal(f1), 
//         //                           tmp_solver.hullGeometry->faceNormal(f2), 
//         //                           Vector3({0,2,0}), 1., 12, 1, 0.1, glm::vec3(1., 0., 1.), 0.5);
        
//         Face ff1 = face_last_face[f1];
//         Face ff2 = face_last_face[f2];
//         // for visuals
//         Eigen::Vector3<Scalar> ff1_normal  = face_normals[ff1], // static
//                                ff2_normal  = face_normals[ff2]; // static
//         Eigen::Vector3<Scalar> adj_f1_normal  = face_normals[f1], // static
//                                adj_f2_normal  = face_normals[f2]; // static
//         // std::cout << " ff1: " << ff1.getIndex() << " ff2: " << ff2.getIndex() << std::endl;
//         // std::cout << " with normals: " << ff1_normal.transpose() << "  --  " << ff2_normal.transpose() << std::endl;
//         for (Vertex v: {e.firstVertex(), e.secondVertex()}){ // two ways to go from a saddle edge
            
//             Eigen::Vector3<Scalar> tmp_normal  = edge_bnd_normal;   // changes in the while loop
//             Vertex current_v = v;                            // changes in the while loop
//             Eigen::Vector3<Scalar> v_normal    = (hull_positions.row(current_v.getIndex()) - G.transpose()).normalized(); // changes in the while loop
            
//             // polyscope::registerPointCloud("current V", std::vector<Vector3>{vec_to_GC_vec3(tmp_normal) +Vector3({0,2,0})});
//             // polyscope::show();

//             // ** DONT CONVERT these to template**
//             double ff1_area_sign = ff1_normal.dot(tmp_normal.cross(v_normal)) >= 0 ? 1. : -1.; // ff1 on rhs of bndN->vN; static
//             double ff2_area_sign = ff2_normal.dot(tmp_normal.cross(v_normal)) >= 0 ? 1. : -1.; // ff2 on rhs of bndN->vN; static
//             double adj_f1_area_sign = adj_f1_normal.dot(tmp_normal.cross(v_normal)) >= 0 ? 1. : -1.; // ff1 on rhs of bndN->vN; static
//             double adj_f2_area_sign = adj_f2_normal.dot(tmp_normal.cross(v_normal)) >= 0 ? 1. : -1.; // ff2 on rhs of bndN->vN; static
            
//             ff1_area_sign *= adj_f1_area_sign == ff1_area_sign ? 1. : -1.;
//             ff2_area_sign *= adj_f2_area_sign == ff2_area_sign ? 1. : -1.;

//             Eigen::Vector3<Scalar> next_normal;                     // changes in the while loop
//             int cnt = 0;
//             while (true) {
//                 cnt++;
//                 // std::cout << "   - at vertex " << current_v.getIndex() << std::endl;
//                 // if (cnt > 50){
//                 //     std::cout << "  - too many iterations\n";
//                 //     break;
//                 // }
//                 if (vertex_is_stabilizable[current_v]){
//                     next_normal = v_normal;
//                     Scalar patch_area_ff1 = triangle_patch_signed_area_on_sphere<Scalar>(ff1_normal, tmp_normal, next_normal);
//                     face_region_area[ff1] += ff1_area_sign * patch_area_ff1;
//                     Scalar patch_area_ff2 = triangle_patch_signed_area_on_sphere<Scalar>(ff2_normal, tmp_normal, next_normal);
//                     face_region_area[ff2] += ff2_area_sign * patch_area_ff2;
//                     // draw_arc_on_sphere_static(vec2vec3(tmp_normal), vec2vec3(next_normal), Vector3({0,2,0}), 1., 12, 0, 0.1, glm::vec3(1., 0., 1.), 0.5);
//                     // polyscope::show();         
//                     // for barycentric energy
//                     face_region_boundary_vertices[ff1].insert(current_v);
//                     face_region_boundary_vertices[ff2].insert(current_v);
//                     break; // Exit while; got to a maximum
//                 }
//                 else{
//                     for (Edge neigh_e: current_v.adjacentEdges()){
//                         // printf(" -e %d-", neigh_e.getIndex());
//                         if (neigh_e != current_e && edge_next_vertex[neigh_e] == current_v){ // neigh_e is a source for this vertex
//                             Face f1 = neigh_e.halfedge().face(),
//                                  f2 = neigh_e.halfedge().twin().face();
//                             bool sign_change = false;
//                             next_normal = intersect_arcs(v_normal, tmp_normal, 
//                                                          face_normals[f1], 
//                                                          face_normals[f2]);
//                             if(next_normal.norm() != 0.){ // found the next source normal
//                                 Scalar patch_area_ff1 = triangle_patch_signed_area_on_sphere<Scalar>(ff1_normal, tmp_normal, next_normal);
//                                 face_region_area[ff1] += ff1_area_sign * patch_area_ff1;
//                                 Scalar patch_area_ff2 = triangle_patch_signed_area_on_sphere<Scalar>(ff2_normal, tmp_normal, next_normal);
//                                 face_region_area[ff2] += ff2_area_sign * patch_area_ff2;

//                                 // //
//                                 // draw_arc_on_sphere_static(vec_to_GC_vec3(tmp_normal), vec_to_GC_vec3(next_normal), Vector3({0,2,0}), 1., 12, 0, 0.1, glm::vec3(1., 0., 1.), 0.5);
//                                 // polyscope::show();
//                                 // updates
//                                 tmp_normal = next_normal;
//                                 current_v = neigh_e.otherVertex(current_v);
//                                 current_e = neigh_e;
//                                 v_normal  = (hull_positions.row(current_v.getIndex()) - G.transpose()).normalized();
//                                 break; // don't check other edges; go for next cell
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }
    

//     //// separate energy terms
//     Scalar total_coplanar_energy = 0.,
//         pure_dice_energy = 0.,
//         total_bary_energy = 0.,
//         total_cluster_distance_energy = 0.;
//     size_t non_zero_clusters = 0;
//     Scalar energy = 0.;
//     if (policy_general == "fair"){
//         // sort for dice energy
//         Scalar total_area = 0.;
//         std::vector<std::pair<Face, Scalar>> probs;
//         for (Face f: tmp_solver.hullMesh->faces()){
//             if (face_region_area[f] > 0){
//                 // face_region_area[f] /= (4.*PI);
//                 probs.push_back({f, face_region_area[f]});
//                 total_area += face_region_area[f];
//             }
//         }
//         // std::cout << "sum over probabilities: " << total_area/(4.*PI) << std::endl;
//         std::sort(probs.begin(), probs.end(), [] (auto a, auto b) { return a.second > b.second; });
//         // printf(" static sorted probs: \n");
//         size_t tmp_side_cnt = 0;
//         double goal_prob = 1. / (double)side_count;
//         for (auto pair: probs){ 
//             Face f = pair.first; // f is stable by construction
//             if (tmp_side_cnt < side_count){ // goal 1/n

//                 Scalar current_cluster_accum_area = 0.;
//                 // pure dice && coplanar energy
//                 for (Face f2: tmp_solver.hullMesh->faces()){
//                     if (tmp_solver.face_last_face[f2] == f2 && // f2 is stable
//                         (face_normals[f2] - face_normals[f]).norm() < 1e-1){ // f2 is close enough to f
//                         // coplanar energy
//                         total_coplanar_energy += (face_normals[f2] - face_normals[f]).squaredNorm();
//                         // dice energy
//                         current_cluster_accum_area += face_region_area[f2];
//                     }
//                 }
//                 // pure dice energy
//                 Scalar diff = current_cluster_accum_area - goal_prob * 4. * PI;
//                 pure_dice_energy += diff * diff;

//                 // barycentric energy term
//                 total_bary_energy += single_DAG_cluster_bary_e<Scalar>(f, face_normals[f], hull_positions, G, face_region_boundary_vertices);
//             }
//             else { // goal zero
//                 Scalar diff = face_region_area[pair.first];
//                 pure_dice_energy += diff * diff;
//             }
//             tmp_side_cnt++;
//         }
//         if (tmp_side_cnt < side_count){
//             // printf("  - not enough faces\n");
//             pure_dice_energy += (side_count - tmp_side_cnt) * goal_prob * goal_prob * 4. * PI * 4. * PI;
//         }
//         energy += pure_dice_energy + bary_reg * total_bary_energy + co_planar_reg * total_coplanar_energy;
//     }
//     else if (policy_general == "manual"){
//         FaceData<double> goal_probs = manual_stable_only_face_prob_assignment(&tmp_solver, normal_prob_assignment);
//         for (Face f: tmp_solver.hullMesh->faces()){
//             // // norm 2
//             Scalar diff = face_region_area[f] - goal_probs[f] * 4. * PI;
//             pure_dice_energy += diff * diff;
//             // // KL divergence
//             //     double goal_area = my_probs[f] * 4. * PI;
//             //     energy += goal_area * log(goal_area / face_region_area[f]);
//             // }

//             // barycentric energy term
//             if (face_last_face[f] == f)
//                 non_zero_clusters++;
//             if (face_region_area[f] > 0 && goal_probs[f] > 0){ // stable face that is supposed to be stable
//                 total_bary_energy += single_DAG_cluster_bary_e<Scalar>(f, face_normals[f], hull_positions, G, face_region_boundary_vertices);
//             }

//             // penalize non-coplanar faces in the rolling DAG
//             if (face_last_face[f] != f){
//                 Eigen::Vector3<Scalar> f_normal = face_normals[f];
//                 Eigen::Vector3<Scalar> last_f_normal = face_normals[face_last_face[f]];
//                 total_coplanar_energy += (f_normal- last_f_normal).squaredNorm();
//             }
//         }
//         energy += pure_dice_energy + co_planar_reg * total_coplanar_energy + bary_reg * total_bary_energy;
//     }
//     else if (policy_general == "manualCluster"){
//         std::vector<std::tuple<std::vector<Face>, double, Vector3>> clustered_probs = 
//                     manual_clustered_face_prob_assignment(&tmp_solver, normal_prob_assignment);
        
//         std::vector<Eigen::Vector3<Scalar>> cluster_barycenters;

//         for (auto cluster: clustered_probs){
//             std::vector<Face> cluster_faces = std::get<0>(cluster);
//             double cluster_prob = std::get<1>(cluster);
//             if (cluster_faces.size() > 0)
//                 non_zero_clusters++;
//             else
//                 continue;

//             // Pure dice energy
//             Scalar cluster_area_sum = 0.;
//             for (Face f: cluster_faces){ // sum of stable face areas in the cluster
//                 if (face_region_area[f] > 0){
//                     cluster_area_sum += face_region_area[f];
//                 }
//             }
//             Scalar diff = cluster_area_sum - cluster_prob * 4. * PI;
//             pure_dice_energy += diff * diff;

//             // Co-planar energy
//             total_coplanar_energy += single_cluster_coplanar_e<Scalar>(cluster_faces, face_normals, face_last_face, unstable_attaction_thresh, verbose);
            
//             // Barycentric energy
//             Eigen::Vector3<Scalar> stable_cluster_faces_barycenter = Eigen::Vector3<Scalar>::Zero(),
//                                    cluster_barycenter = Eigen::Vector3<Scalar>::Zero();
//             std::tie(cluster_barycenter, stable_cluster_faces_barycenter) = 
//             region_and_stables_barycenters(cluster_faces, face_last_face, 
//                                            hull_positions, G,
//                                            face_normals, face_region_boundary_vertices);
//             total_bary_energy += (cluster_barycenter - stable_cluster_faces_barycenter).squaredNorm();
            
//             // cluster distance term
//             cluster_barycenters.push_back(stable_cluster_faces_barycenter);
//         }
//         // cluster distance term
//         for (int i = 0; i < cluster_barycenters.size(); i++){
//             for (int j = i+1; j < cluster_barycenters.size(); j++){
//                 Scalar diff_sqr = (cluster_barycenters[i] - cluster_barycenters[j]).squaredNorm();
//                 total_cluster_distance_energy += 1./(diff_sqr);
//             }
//         }
//         // total energy
//         energy += pure_dice_energy + 
//                   co_planar_reg * total_coplanar_energy + bary_reg * total_bary_energy + cluster_distance_reg * total_cluster_distance_energy;
//     }
//     else {
//         throw std::invalid_argument("policy not recognized");
//     }
//     if (verbose){
//         if (bary_reg > 0)             std::cout << "  - barycenter dist energy : " << bary_reg << " * " <<  total_bary_energy << std::endl;
//         if (co_planar_reg > 0)        std::cout << "  - coplanar energy        : " << co_planar_reg << " * " << total_coplanar_energy << std::endl;
//         if (cluster_distance_reg > 0) std::cout << "  - cluster distance energy: " << cluster_distance_reg << " * " << total_cluster_distance_energy << std::endl;
//         std::cout << "  - pure dice energy: " << pure_dice_energy << "\t Nz cluster count: " << non_zero_clusters << std::endl;
//         std::cout << "\t\t\t** total energy: " << energy << std::endl;
//     }
//     return energy;
// }

// template <typename Scalar>
// Eigen::Vector3<Scalar> BoundaryBuilder::point_to_segment_normal(Eigen::Vector3<Scalar> P, Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B){
//     Eigen::Vector3<Scalar> PB = B - P,
//                     AB = B - A;
//     return (PB * AB.dot(AB) - AB * AB.dot(PB)).normalized();
// }

// template <typename Scalar>
// Eigen::Vector3<Scalar> BoundaryBuilder::intersect_arcs(Eigen::Vector3<Scalar> v_normal, Eigen::Vector3<Scalar> R2, Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B){
//     Eigen::Vector3<Scalar> p = (((v_normal.cross(R2)).cross(A.cross(B)))).normalized();
//     if (p.cross(B).dot(A.cross(B)) >= 0. && p.cross(A).dot(B.cross(A)) >= 0.)
//         return p;
//     else if ((-p).cross(B).dot(A.cross(B)) >= 0. && (-p).cross(A).dot(B.cross(A)) >= 0.)
//         return -p;
//     else
//         return Eigen::Vector3<Scalar>::Zero();
//     // if (p.dot(A) >= A.dot(B) && p.dot(B) >= A.dot(B))
//     //     return p;
//     // else if (-p.dot(A) >= A.dot(B) && -p.dot(B) >= A.dot(B))
//     //     return -p;
//     // else 
//     //     return Eigen::Vector3<Scalar>::Zero();
// }

// template <typename Scalar>
// Scalar BoundaryBuilder::triangle_patch_signed_area_on_sphere(Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B, Eigen::Vector3<Scalar> C){
//     Scalar Adotbc = A.dot(B.cross(C));
//     Scalar denom = A.norm()*B.norm()*C.norm() + A.dot(B)*C.norm() + B.dot(C)*A.norm() + C.dot(A)*B.norm();
//     return 2. * atan2(Adotbc, denom);
// }


// template <typename Scalar>
// Scalar BoundaryBuilder::single_cluster_coplanar_e(std::vector<Face> faces, 
//                                         FaceData<Eigen::Vector3<Scalar>> face_normals,
//                                         FaceData<Face> face_last_face, double unstable_attaction_thresh,
//                                         bool verbose){
//     Scalar local_energy = 0.;
//     for (Face f1: faces){
//         Face ff1 = face_last_face[f1];
//         if (ff1 == f1){ // Could remove
//             if (std::find(faces.begin(), faces.end(), ff1) != faces.end()){ // final face is in the same cluster
//                 for (Face f2: faces){
//                     Face ff2 = face_last_face[f2];
//                     if (ff2 == f2 && f1 != f2){ // Could remove
//                         if (std::find(faces.begin(), faces.end(), ff2) != faces.end()){ // final face is in the same cluster
//                             local_energy += (face_normals[f1] - face_normals[f2]).squaredNorm();
//                         }
//                     }
//                 }
//             }
//         }
//         else {
//             if ((face_normals[f1] - face_normals[ff1]).norm() < unstable_attaction_thresh){ // close enough unstab // 0.1 up to 6 prism
//                 // if (verbose)
//                 //     std::cout << "  - f" << f1.getIndex() << "(" << f1.halfedge().vertex().getIndex() << ", " << f1.halfedge().tipVertex().getIndex() << ", " << f1.halfedge().next().tipVertex().getIndex() << ")" 
//                 //     << " is close to f" << ff1.getIndex()<< "("<< ff1.halfedge().vertex().getIndex() << ", " << ff1.halfedge().tipVertex().getIndex() << ", " << ff1.halfedge().next().tipVertex().getIndex() << ")"
//                 //     << std::endl;
//                 local_energy += (face_normals[f1] - face_normals[ff1]).squaredNorm();
//             }
//         }
//     }
//     // std::vector<std::pair<size_t, size_t>> face_ind_pairs{{0,4}, {2,3}, {7,8}};
//     // for (auto pair: face_ind_pairs){
//     //     size_t f1 = pair.first,
//     //            f2 = pair.second;
//     //     local_energy += (face_normals[f1] - face_normals[f2]).squaredNorm();
//     // }
//     return local_energy;
// }


// template <typename Scalar> 
// Scalar BoundaryBuilder::single_cluster_bary_e(std::vector<Face> faces, FaceData<Face> face_last_face,
//                                             Eigen::MatrixX3<Scalar> hull_positions, Eigen::Vector3<Scalar> G,
//                                             FaceData<Eigen::Vector3<Scalar>> face_normals,
//                                             FaceData<std::set<Vertex>> face_region_boundary_vertices){
//     Eigen::Vector3<Scalar> stable_cluster_faces_barycenter = Eigen::Vector3<Scalar>::Zero();
//     Eigen::Vector3<Scalar> cluster_barycenter = Eigen::Vector3<Scalar>::Zero();

//     std::tie(cluster_barycenter, stable_cluster_faces_barycenter) = 
//     region_and_stables_barycenters(faces, face_last_face, 
//                                    hull_positions, G, 
//                                    face_normals, face_region_boundary_vertices);
//     Scalar stable_face_distances_to_barycenter = 0.;
    
//     // -- Penalize very stable face distance to bary center; could de-motivate splitting; // sum (fi - bary)^2
//     // for (Face cluster_f: cluster_faces)
//     //     if (face_region_area[cluster_f] > 0) // stable cluster faces
//     //         stable_face_distances_to_barycenter += (cluster_barycenter - face_normals[cluster_f]).squaredNorm();
    
//     // -- Penalize only stable barycenters distance to cluster barycenter; might leave splits alone // [mean (fi) - bary]^2
//     stable_face_distances_to_barycenter = (cluster_barycenter - stable_cluster_faces_barycenter).squaredNorm();
    
//     return stable_face_distances_to_barycenter;
// }


// template <typename Scalar>
// std::pair<Eigen::Vector3<Scalar>, Eigen::Vector3<Scalar>> 
// BoundaryBuilder::region_and_stables_barycenters(std::vector<Face> faces, FaceData<Face> face_last_face,
//                               Eigen::MatrixX3<Scalar> hull_positions, Eigen::Vector3<Scalar> G,
//                               FaceData<Eigen::Vector3<Scalar>> face_normals,
//                               FaceData<std::set<Vertex>> face_region_boundary_vertices){
//     Eigen::Vector3<Scalar> stable_cluster_faces_barycenter = Eigen::Vector3<Scalar>::Zero();
//     Eigen::Vector3<Scalar> cluster_barycenter = Eigen::Vector3<Scalar>::Zero();
//     std::set<Vertex> cluster_boundary_vertices;
    
//     double stable_face_count = 0;
//     if (faces.size() > 0){
//         for (Face cluster_f: faces){
//             if (face_last_face[cluster_f] == cluster_f){ // stable face in cluster
//                 cluster_boundary_vertices.insert(face_region_boundary_vertices[cluster_f].begin(), face_region_boundary_vertices[cluster_f].end());
//                 stable_cluster_faces_barycenter += face_normals[cluster_f];
//                 stable_face_count += 1.;
//             }
//         }
//         if (stable_face_count > 0){
//             stable_cluster_faces_barycenter = stable_cluster_faces_barycenter.normalized();
//             for (Vertex v: cluster_boundary_vertices)
//                 cluster_barycenter += (hull_positions.row(v.getIndex()) - G.transpose()).normalized(); // vertex-G normals
//             cluster_barycenter /= (double)cluster_boundary_vertices.size();
//             cluster_barycenter = cluster_barycenter.normalized();
//         }
//     }

//     return {cluster_barycenter, stable_cluster_faces_barycenter};
// }


// // --- for Fair and Manual policies ---
// template <typename Scalar> 
// Scalar 
// BoundaryBuilder::single_DAG_cluster_bary_e(Face stable_face, Eigen::Vector3<Scalar> stable_face_normal,
//                         Eigen::MatrixX3<Scalar> hull_positions, Eigen::Vector3<Scalar> G,
//                         FaceData<std::set<Vertex>> face_region_boundary_vertices){
//     Eigen::Vector3<Scalar> barycenter = Eigen::Vector3<Scalar>::Zero();
//     for (Vertex v: face_region_boundary_vertices[stable_face]){
//         barycenter += (hull_positions.row(v.getIndex())- G.transpose()).normalized(); // vertex-G normals
//     }
//     barycenter /= face_region_boundary_vertices[stable_face].size();
//     barycenter = barycenter.normalized();
//     return (barycenter - stable_face_normal).squaredNorm(); // diff^2
// }

// #endif