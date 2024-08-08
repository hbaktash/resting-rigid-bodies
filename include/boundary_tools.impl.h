// Forward declarations
void draw_arc_on_sphere_static(Vector3 p1, Vector3 p2, Vector3 center, double radius, 
                        size_t seg_count, size_t edge_ind, 
                        double radi_scale, glm::vec3 color, 
                        float arc_curve_radi);

template <typename Scalar>
Scalar BoundaryBuilder::dice_energy(Eigen::MatrixX3<Scalar> hull_positions, Eigen::Vector3<Scalar> G,
                                    Forward3DSolver &tmp_solver, size_t side_count
                                    ){
    // precomputes
    // printf("initalize precomputes\n");
    tmp_solver.updated = false;
    tmp_solver.initialize_pre_computes();
    FaceData<Face> face_last_face = tmp_solver.face_last_face;
    VertexData<bool> vertex_is_stabilizable = tmp_solver.vertex_is_stabilizable; 
    EdgeData<Vertex> edge_next_vertex = tmp_solver.edge_next_vertex;
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
    // printf("  back-flowing terminal edges %lu \n", terminal_edges.size());
    int i = 1;
    for (Edge e: terminal_edges){
        // printf("\n - starting at terminal edge: %d/%d \n", i++, terminal_edges.size());
        // std::cout << "  -e" << e.getIndex() << " : " << e.firstVertex().getIndex() << " , " << e.secondVertex().getIndex() << " adj faces: "<< e.halfedge().face().getIndex() << " " << e.halfedge().twin().face().getIndex() << std::endl;
        // bnd_normal->normal_ad = point_to_segment_normal_ad(e);
        Eigen::Vector3<Scalar> edge_bnd_normal = point_to_segment_normal<Scalar>(G, hull_positions.row(e.firstVertex().getIndex()), hull_positions.row(e.secondVertex().getIndex()));
        // std::vector<std::pair<size_t, size_t>> vertex_pairs;
        // vertex_pairs.push_back({e.firstVertex().getIndex(), e.secondVertex().getIndex()});
        // polyscope::registerCurveNetwork("current edge", hull_positions, vertex_pairs);
        Edge current_e = e;
        Face f1 = current_e.halfedge().face(),
             f2 = current_e.halfedge().twin().face();
        // Eigen::Vector3<Scalar> e0{1,0,0};
        // double f10 = face_normals[f1].dot(),
        //        f11 = face_normals[f1].dot(Eigen::Vector3<Scalar>{0,1,0}),
        //        f12 = face_normals[f1].dot(Eigen::Vector3<Scalar>{0,1,0});
        // auto f11 = face_normals[f1].coef(0);//, f12 = face_normals[f1](2), f20 = face_normals[f2](0), f21 = face_normals[f2](1), f22 = face_normals[f2](2);
        // draw_arc_on_sphere_static(Vector3({f10, f11, f12}), Vector3({f20, f21, f22}), Vector3({0,2,0}), 1., 12, 1, 0.1, glm::vec3(1., 0., 1.), 0.5);
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
                //     // std::cout << "  - too many iterations\n";
                //     break;
                // }
                if (vertex_is_stabilizable[current_v]){
                    next_normal = v_normal;
                    Scalar patch_area_ff1 = triangle_patch_signed_area_on_sphere<Scalar>(ff1_normal, tmp_normal, next_normal);
                    face_region_area[ff1] += ff1_area_sign * patch_area_ff1;
                    Scalar patch_area_ff2 = triangle_patch_signed_area_on_sphere<Scalar>(ff2_normal, tmp_normal, next_normal);
                    face_region_area[ff2] += ff2_area_sign * patch_area_ff2;
                    // draw_arc_on_sphere_static(Vector3({tmp_normal(0),tmp_normal(1),tmp_normal(2)}), Vector3({next_normal(0),next_normal(1),next_normal(2)}), Vector3({0,2,0}), 1., 12, 0, 0.1, glm::vec3(1., 0., 1.), 0.5);
                    // polyscope::show();            
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
    Scalar energy = 0.;
    size_t tmp_side_cnt = 0;
    double goal_prob = 1. / (double)side_count;
    for (auto pair: probs){ 
        // std::cout << "  -f" << pair.first.getIndex() << " : " << pair.second/(4. * PI) << std::endl; 
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
    
    return energy; 
}

template <typename Scalar>
Eigen::Vector3<Scalar> BoundaryBuilder::point_to_segment_normal(Eigen::Vector3<Scalar> P, Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B){
    Eigen::Vector3<Scalar> PB = B - P,
                    AB = B - A;
    return (PB * AB.dot(AB) - AB * AB.dot(PB)).normalized();
}
// autodiff::Vector3var PB = var_positions.row(v2.getIndex()) - var_G.transpose(),// var_G,
//                          AB = var_positions.row(v2.getIndex()) - var_positions.row(v1.getIndex());
//     return (PB * AB.dot(AB) - AB * AB.dot(PB));//.normalized();

template <typename Scalar>
Eigen::Vector3<Scalar> BoundaryBuilder::intersect_arcs(Eigen::Vector3<Scalar> v_normal, Eigen::Vector3<Scalar> R2, Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B){
    Eigen::Vector3<Scalar> p = (((v_normal.cross(R2)).cross(A.cross(B)))).normalized(); // var_G.transpose()
    if (p.dot(A) >= A.dot(B) && p.dot(B) >= A.dot(B))
        return p;
    else if (-p.dot(A) >= A.dot(B) && -p.dot(B) >= A.dot(B))
        return -p;
    else 
        return Eigen::Vector3<Scalar>::Zero();
}

template <typename Scalar>
Scalar BoundaryBuilder::triangle_patch_signed_area_on_sphere(Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B, Eigen::Vector3<Scalar> C){
    Scalar Adotbc = A.dot(B.cross(C));
    Scalar denom = A.norm()*B.norm()*C.norm() + A.dot(B)*C.norm() + B.dot(C)*A.norm() + C.dot(A)*B.norm();
    return 2. * atan2(Adotbc, denom);
}