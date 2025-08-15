#include "boundary_tools.h"


// // temp visuals
// void tmp_arc_vis(Vector3 p1, Vector3 p2, Vector3 center, double radius, 
//                  size_t seg_count, size_t edge_ind, 
//                  double radi_scale, glm::vec3 color){
//   // p1, p2 just represent normal vectors
//   if (norm(p1) > 1.01)
//     polyscope::warning("p1 norm larger than 1!");
//   if (norm(p2) > 1.01)
//     polyscope::warning("p2 norm larger than 1!");

//   std::vector<std::array<size_t, 2>> edgeInds;
//   std::vector<Vector3> positions;
//   double sqrt_radi = sqrt(radius);
//   // walk on p1-p2 segment
//   Vector3 curr_point = p1,
//           forward_vec = (p2-p1)/(double)seg_count;
//   Vector3 next_point = curr_point + forward_vec;
//   Vector3 curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center ,
//           next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
//   positions.push_back(curr_point_on_sphere);
//   for (size_t i = 0; i < seg_count; i++){
//     // add to positions list
//     curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center ,
//     next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
//     positions.push_back(next_point_on_sphere);
//     // add segment indices
//     edgeInds.push_back({i, i+1});

//     // update points
//     curr_point = next_point;
//     next_point += forward_vec;
//   }
//   auto psArcCurveNet = polyscope::registerCurveNetwork("Arc curve " + std::to_string(edge_ind), positions, edgeInds);
//   psArcCurveNet->setRadius(0.1 * radi_scale, false);
//   psArcCurveNet->setColor(color);
//   psArcCurveNet->setEnabled(true);
// }

size_t BoundaryNormal::counter = 0;
// constructor
BoundaryNormal::BoundaryNormal(Vector3 _normal)
    : index(counter++){
        normal = _normal;
        host_e = Edge();
        host_v = Vertex();
}

void BoundaryNormal::add_neighbor(BoundaryNormal* _neigh) {
    neighbors.push_back(_neigh);
}

// constructor
BoundaryBuilder::BoundaryBuilder(Forward3DSolver *forward_solver_){
    forward_solver = forward_solver_;
    // mesh = forward_solver->hullMesh;
    // geometry = forward_solver->hullGeometry;
}

std::vector<Edge> BoundaryBuilder::find_terminal_edges(){
    std::vector<Edge> terminal_edges;
    // printf("  buidling face-last-face\n");
    // printf("  finding terminal edges \n");
    size_t cnt = 0;
    for (Edge e: forward_solver->hullMesh->edges()){
        if (forward_solver->edge_next_vertex[e].getIndex() == INVALID_IND){ // singular edge
            Face f1 = e.halfedge().face(),
                 f2 = e.halfedge().twin().face();
            if (forward_solver->face_last_face[f1] != forward_solver->face_last_face[f2]){ // saddle edge
                // proved: at least one singular edge like this must exist
                // TODO: assert that stable normal falls inside the edge arc 
                terminal_edges.push_back(e);
                cnt++;
            }
        }
    }
    return terminal_edges;
}

void BoundaryBuilder::build_boundary_normals(){
    vertex_boundary_normal = VertexData<BoundaryNormal*>(*forward_solver->hullMesh, nullptr);
    edge_boundary_normals  = EdgeData<std::vector<BoundaryNormal*>>(*forward_solver->hullMesh);
    // edge_boundary_normals = EdgeData<std::vector<Vector3>>(*forward_solver->hullMesh);
    BoundaryNormal::counter = 0;
    // 
    // printf("  buidling face-last-face\n");
    // printf("  finding terminal edges \n");
    std::vector<Edge> terminal_edges = find_terminal_edges();

    // final ascending manifold areas
    face_region_area = FaceData<double>(*forward_solver->hullMesh, 0.);
    
    // 
    MS_neighbor_distances = FaceData<FaceData<double>>(*forward_solver->hullMesh);
    for (Face f: forward_solver->hullMesh->faces()){
        MS_neighbor_distances[f] = FaceData<double>(*forward_solver->hullMesh, -1.);
    }
    // back-flow from all terminal edges
    // printf("  back-flowing terminal edges \n");
    int i = 0;
    for (Edge e: terminal_edges){
        Vector3 stable_edge_normal = forward_solver->edge_stable_normal[e];
        BoundaryNormal *bnd_normal = new BoundaryNormal(stable_edge_normal);
        Face f1 = e.halfedge().face(),
             f2 = e.halfedge().twin().face();
        bnd_normal->f1 = forward_solver->face_last_face[f1];
        bnd_normal->f2 = forward_solver->face_last_face[f2];
        bnd_normal->host_e = e;
        
        MS_neighbor_distances[bnd_normal->f1][bnd_normal->f2] = angle(forward_solver->hullGeometry->faceNormal(bnd_normal->f1), 
                                                                      stable_edge_normal);
        MS_neighbor_distances[bnd_normal->f2][bnd_normal->f1] = angle(forward_solver->hullGeometry->faceNormal(bnd_normal->f2), 
                                                                      stable_edge_normal);
        // // DEBUG
        // tmp_arc_vis(forward_solver->hullGeometry->faceNormal(f1), 
        //                           forward_solver->hullGeometry->faceNormal(f2), 
        //                           Vector3({0,2,0}), 1., 12, 1, 0.1, glm::vec3(1., 0., 1.));
        // for visuals
        edge_boundary_normals[e].push_back(bnd_normal);
        
        Vector3 adj_f1_normal  = forward_solver->hullGeometry->faceNormal(f1),
                adj_f2_normal  = forward_solver->hullGeometry->faceNormal(f2);
        for (Vertex v: {e.firstVertex(), e.secondVertex()}){
            Vector3 tmp_normal = bnd_normal->normal,
                    ff1_normal  = forward_solver->hullGeometry->faceNormal(bnd_normal->f1), // final faces
                    ff2_normal  = forward_solver->hullGeometry->faceNormal(bnd_normal->f2),
                    v_normal   = forward_solver->vertex_stable_normal[v];
            // Vector3 imm_f1_normal = forward_solver->hullGeometry->faceNormal(e.halfedge().face()), // immediate face neighbors
            //         imm_f2_normal = forward_solver->hullGeometry->faceNormal(e.halfedge().twin().face());
            
            // std::cout << " ff1: " << bnd_normal->f1.getIndex() << " ff2: " << bnd_normal->f2.getIndex() << std::endl;
            // std::cout << " with normals: " << ff1_normal.transpose() << "  --  " << ff2_normal.transpose() << std::endl;

            // OLD
            // double f1_area_sign = dot(ff1_normal, cross(v_normal, tmp_normal)) >= 0 ? 1. : -1.; // f1 on rhs of bndN->vN
            
            double ff1_area_sign = dot(ff1_normal, cross(tmp_normal,v_normal)) >= 0 ? 1. : -1.; // ff1 on rhs of bndN->vN; static
            double ff2_area_sign = dot(ff2_normal, cross(tmp_normal,v_normal)) >= 0 ? 1. : -1.; // ff2 on rhs of bndN->vN; static
            double adj_f1_area_sign = dot(adj_f1_normal, cross(tmp_normal, v_normal)) >= 0 ? 1. : -1.; // ff1 on rhs of bndN->vN; static
            double adj_f2_area_sign = dot(adj_f2_normal, cross(tmp_normal, v_normal)) >= 0 ? 1. : -1.; // ff2 on rhs of bndN->vN; static
            
            ff1_area_sign *= adj_f1_area_sign == ff1_area_sign ? 1. : -1.;
            ff2_area_sign *= adj_f2_area_sign == ff2_area_sign ? 1. : -1.;

            if (forward_solver->vertex_is_stabilizable[v])
                flow_back_boundary_on_edge(bnd_normal, Edge(), v, 
                                           ff1_area_sign, ff2_area_sign);
            else {
                for (Edge neigh_e: v.adjacentEdges()){
                    if (neigh_e != e &&
                        forward_solver->edge_next_vertex[neigh_e] == v){ // neigh_e is a source for this vertex
                        // printf("-SE %d-", neigh_e.getIndex());
                        // face_attraction_boundary[bnd_normal->f1].push_back(bnd_normal);
                        // face_attraction_boundary[bnd_normal->f2].push_back(bnd_normal);
                        flow_back_boundary_on_edge(bnd_normal, neigh_e, v, 
                                                ff1_area_sign, ff2_area_sign);
                    }
                }
            }
            // printf(" \n One side done! \n");
        }
    }
}


// recursively follow the boundary curve to a source
bool BoundaryBuilder::flow_back_boundary_on_edge(BoundaryNormal* bnd_normal, Edge src_e, Vertex common_vertex,
                                                 double f1_area_sign, double f2_area_sign){
    // bnd_normal has to be in boundary normals of dest_e; won't assert tho for better performance
    Vertex v = common_vertex; // given as argument for better performance
    
    
    Vector3 tmp_normal({0.,0.,0.}), next_normal({0.,0.,0.});
    if (forward_solver->vertex_is_stabilizable[v]){ // we are at a source
        BoundaryNormal *curr_vertex_boundary_normal = vertex_boundary_normal[v];
        if (vertex_boundary_normal[v] == nullptr){ // first time arriving at this stable vertex; create the boundary normal
            Vector3 stable_vertex_normal = forward_solver->vertex_stable_normal[v];
            curr_vertex_boundary_normal = new BoundaryNormal(stable_vertex_normal);
            vertex_boundary_normal[v] = curr_vertex_boundary_normal;
            curr_vertex_boundary_normal->host_v = v;
        }
        curr_vertex_boundary_normal->add_neighbor(bnd_normal);
        bnd_normal->add_neighbor(curr_vertex_boundary_normal);
        // handling face region boundary
        // face_attraction_boundary[bnd_normal->f1].push_back(curr_vertex_boundary_normal);
        // face_attraction_boundary[bnd_normal->f2].push_back(curr_vertex_boundary_normal);

        tmp_normal = bnd_normal->normal;
        next_normal = curr_vertex_boundary_normal->normal;
        // printf("  -> V done! \n ");

        // // DEBUG
        // std::cout << "at stable V: " << v.getIndex() << std::endl;
        // polyscope::registerPointCloud("Stable V", std::vector<Vector3>{next_normal +Vector3({0,2,0})})->setPointRadius(0.02);
        // tmp_arc_vis(tmp_normal, next_normal, Vector3({0,2,0}), 1., 12, 0, 0.1, glm::vec3(1., 0., 1.));
        // polyscope::show();  
    }
    else {
        // vertex is not an equilibria; i.e. source is outside the vertex
        if (forward_solver->edge_next_vertex[src_e] == v){ // src_e is a source for this vertex
            Face f1 = src_e.halfedge().face(),
                 f2 = src_e.halfedge().twin().face();
                
            // ** actually the following condition does not necessarily have to hold
            //    commented:
            // if (forward_solver->face_last_face[f1] == forward_solver->face_last_face[f2])
            //     return; // this source edge is not a source for this boundary normal 

            Vector3 vertex_stable_normal = forward_solver->vertex_stable_normal[v],
                    tmp_f1_normal = forward_solver->hullGeometry->faceNormal(f1),
                    tmp_f2_normal = forward_solver->hullGeometry->faceNormal(f2);
            bool sign_change = false;
            Vector3 e_bnd_normal = intersect_arc_ray_with_arc(vertex_stable_normal, bnd_normal->normal, 
                                                              tmp_f1_normal, tmp_f2_normal, sign_change);
            // printf("at tmp_e %d, %d. side v %d. e: %d, %d \n", dest_e.firstVertex().getIndex(), dest_e.secondVertex().getIndex(), v.getIndex(), e.firstVertex().getIndex(), e.secondVertex().getIndex());
            if(e_bnd_normal.norm() == 0.){
                // printf(" -src but miss- ");
                return false; // not a source for the given bnd_normal
            }
            // found the source normal
            // printf(" -> e %d ", src_e.getIndex());

            BoundaryNormal *new_boundary_normal = new BoundaryNormal(e_bnd_normal);
            edge_boundary_normals[src_e].push_back(new_boundary_normal);
            // edge_boundary_normals[src_e].push_back(e_bnd_normal);

            bnd_normal->add_neighbor(new_boundary_normal);
            new_boundary_normal->add_neighbor(bnd_normal);
            new_boundary_normal->host_e = src_e;

            // region labels
            new_boundary_normal->f1 = bnd_normal->f1;
            new_boundary_normal->f2 = bnd_normal->f2;

            tmp_normal  = bnd_normal->normal,
            next_normal = new_boundary_normal->normal;
            
            // std::cout << "intersected these to get: \n";
            // tmp_arc_vis(vertex_stable_normal, bnd_normal->normal, Vector3({0,2,0}), 1., 12, 0, 0.1, glm::vec3(.1, 1., 0.1));
            // tmp_arc_vis(tmp_f1_normal, tmp_f2_normal, Vector3({0,2,0}), 1., 12, 1, 0.1, glm::vec3(.1, .1, 1.));
            // polyscope::registerPointCloud("4 intersections", std::vector<Vector3>{vertex_stable_normal +Vector3({0,2,0}),
            //                                                                       bnd_normal->normal   +Vector3({0,2,0}),
            //                                                                       tmp_f1_normal        +Vector3({0,2,0}), 
            //                                                                       tmp_f2_normal        +Vector3({0,2,0})})->setPointRadius(0.01);
            // polyscope::show();
            // std::cout << "going through current vertex: " << v.getIndex() << std::endl;
            // polyscope::registerPointCloud("current V", std::vector<Vector3>{vertex_stable_normal +Vector3({0,2,0})})->setPointRadius(0.01);
            // printf("at next bnd normal : \n");//DEBUG
            // tmp_arc_vis(tmp_normal, next_normal, Vector3({0,2,0}), 1., 12, 0, 0.1, glm::vec3(1., 0., 1.));
            // polyscope::show();

            // go with the back-flow
            Vertex next_v = src_e.otherVertex(v);
            if (forward_solver->vertex_is_stabilizable[next_v]){
                // printf(" v ");
                flow_back_boundary_on_edge(new_boundary_normal, Edge(), next_v, f1_area_sign, f2_area_sign);
            }
            else {
                for (Edge next_src_e: next_v.adjacentEdges()){
                    if (next_src_e != src_e &&
                        forward_solver->edge_next_vertex[next_src_e] == next_v){ // next_src_e is a source for this next vertex
                        // std::cout << " e " << next_src_e.getIndex() << std::endl;
                        // non_singularity and divisive-ness will be checked inside the function
                                                
                        bool res = flow_back_boundary_on_edge(new_boundary_normal, next_src_e, next_v, 
                                                              f1_area_sign, f2_area_sign);
                        if (res) // got to maximum and returning
                            break;
                    }
                }
            }
        }
        else{
            // std::cout << " not src! " << std::endl;
            ;
        }
    }
    // std::cout << "here???\n" << std::endl;
    if (tmp_normal.norm() == 0.){ // when does this happen?
        // std::cout << "ret\n" << std::endl;
        return false;
    }
    Vector3 f1_normal = forward_solver->hullGeometry->faceNormal(bnd_normal->f1), // f1,f2 are the same along the current path to maximum
            f2_normal = forward_solver->hullGeometry->faceNormal(bnd_normal->f2);
    // double curr_f1_alignment = dot(f1_normal, cross(next_normal, tmp_normal)) >= 0 ? 1. : -1.; // checking alignment again since it could change along the way
    // double curr_f2_alignment = dot(f2_normal, cross(next_normal, tmp_normal)) >= 0 ? 1. : -1.;
    // double f1_sign_change = f1_area_sign == curr_f1_alignment ? 1. : -1;
    // double f2_sign_change = (-f1_area_sign == curr_f2_alignment) ? 1.: -1;
    face_region_area[bnd_normal->f1] += f1_area_sign * triangle_patch_area_on_sphere(f1_normal, bnd_normal->normal, next_normal);
    face_region_area[bnd_normal->f2] += f2_area_sign * triangle_patch_area_on_sphere(f2_normal, bnd_normal->normal, next_normal);
    // printf(" -ret- ");
    return true;
}


double BoundaryBuilder::get_fair_dice_energy(size_t side_count){
    double energy = 0., goal_area = 4.*PI/(double)side_count, goal_prob = 1./(double)side_count;
    std::vector<double> face_areas;
    for (Face f: forward_solver->hullMesh->faces()){
        if (forward_solver->face_last_face[f] == f)
            face_areas.push_back(face_region_area[f]);
    }
    std::sort(face_areas.begin(), face_areas.end());
    double nth_largest_area = face_areas[face_areas.size() - std::min(side_count, face_areas.size())];
    // Squared norm
    for (Face f: forward_solver->hullMesh->faces()){
        if (forward_solver->face_last_face[f] == f && 
                    face_region_area[f] >= nth_largest_area){
            double diff = goal_area - face_region_area[f];
            energy += diff * diff;
        }
        else if (forward_solver->face_last_face[f] == f){ // noisy stable
            energy += face_region_area[f] * face_region_area[f]; // 0 goal area
        }
    }
    if (side_count > face_areas.size()){
        // printf(" adding penalty for nnot enough stable faces\n");
        energy += ((double)(side_count - face_areas.size())) * goal_area * goal_area; // adding penalty for lack of stable faces
    }
    // KL divergence
    // for (Face f: forward_solver->hullMesh->faces()){
    //     if (forward_solver->face_last_face[f] == f && 
    //                 face_region_area[f] >= nth_largest_area){
    //         double prob = face_region_area[f]/(4.*PI);
    //         energy += prob * std::log(prob/goal_prob);
    //     }
    // }
    return energy;
}

void BoundaryBuilder::print_area_of_boundary_loops(){
    printf(" Face probs:\n");
    std::vector<std::pair<Face, double>> probs;
    double sum = 0.;
    for (Face f: forward_solver->hullMesh->faces()){
        if (face_region_area[f] != 0){
            // face_region_area[f] /= (4.*PI);
            probs.push_back({f, face_region_area[f]/(4.*PI)});
            sum += face_region_area[f]/(4.*PI);
        }
    }
    std::sort(probs.begin(), probs.end(), [] (auto a, auto b) { return a.second > b.second; });
    printf("sorted probs: \n");
    for (auto pair: probs)
        printf("  -f %d (%d,%d,%d): %f\n", pair.first.getIndex(), 
                                           pair.first.halfedge().vertex().getIndex(), 
                                           pair.first.halfedge().tipVertex().getIndex(),
                                           pair.first.halfedge().next().tipVertex().getIndex(),
                                           pair.second);
    printf("sum: %f\n", sum);
}


double BoundaryBuilder::print_pairwise_distances(bool verbose){
  std::vector<std::pair<std::pair<Face, Face>, double>> ratio_pairs;
  for (Face f1: forward_solver->hullMesh->faces()){
    for (Face f2: forward_solver->hullMesh->faces()){
      if (MS_neighbor_distances[f1][f2] != -1.){
        double distance_sum = MS_neighbor_distances[f1][f2] + MS_neighbor_distances[f2][f1];
        // std::cout << "distance sum: " << distance_sum << std::endl;
        if (distance_sum > 0)
          ratio_pairs.push_back({{f1, f2}, std::max({MS_neighbor_distances[f1][f2]/distance_sum, MS_neighbor_distances[f2][f1]/distance_sum})});
      }
    }
  }
  std::sort(ratio_pairs.begin(), ratio_pairs.end(), 
            [&] (auto a, auto b) { return a.second > b.second; });
  if (verbose)
    printf("sorted ratios: \n");
  double area_weighted_ratio_test = 0.;
  for (auto pair: ratio_pairs){
    Face f1 = pair.first.first, 
         f2 = pair.first.second;
    if (verbose){
      std:: cout << "faces: \t" << f1.getIndex() << " \t" << f2.getIndex() << " \t ratio: " << pair.second << 
      " \n\t\tdist:" << MS_neighbor_distances[f1][f2] <<  " \trevert: " << MS_neighbor_distances[f2][f1] << "\n";
    }
    area_weighted_ratio_test += (pair.second - 0.5) * (std::min({face_region_area[f1], face_region_area[f2]}));
  }
  return area_weighted_ratio_test;
}


std::pair<std::vector<std::pair<size_t, size_t>>,std::vector<Vector3>> 
                BoundaryBuilder::MS_complex_edges(){
  std::vector<Vector3> boundary_normals(BoundaryNormal::counter);
  std::set<std::pair<size_t, size_t>> drawn_pairs;
  // printf("showing boundary patches\n");
  // printf("  building pairs\n ");
  for (Edge e: this->forward_solver->hullMesh->edges()){
    for (BoundaryNormal *tmp_bnd_normal: this->edge_boundary_normals[e]){
      if (tmp_bnd_normal != nullptr){
        Vector3 tmp_normal = tmp_bnd_normal->normal;
        boundary_normals[tmp_bnd_normal->index] = tmp_normal; // + gm_shift;
        for (BoundaryNormal *neigh_bnd_normal: tmp_bnd_normal->neighbors){
          Vector3 neigh_normal = neigh_bnd_normal->normal;
          boundary_normals[neigh_bnd_normal->index] = neigh_normal;
          std::pair<size_t, size_t> tmp_pair = {tmp_bnd_normal->index, neigh_bnd_normal->index};
          if (drawn_pairs.find(tmp_pair) == drawn_pairs.end()){
            drawn_pairs.insert(tmp_pair);
          }
        }
      }
    }
  }
  // set to vector
  std::vector<std::pair<size_t, size_t>> ind_pairs_vector;
  ind_pairs_vector.assign(drawn_pairs.begin(), drawn_pairs.end());
  return {ind_pairs_vector, boundary_normals};
}


std::pair<std::vector<std::pair<size_t, size_t>>,std::vector<Vector3>> 
    BoundaryBuilder::MS_complex_edges_of_face(Face f){
  std::vector<Vector3> boundary_normals(BoundaryNormal::counter);
  std::set<std::pair<size_t, size_t>> edge_set;
  // printf("showing boundary patches\n");
  // printf("  building pairs\n ");
  for (Edge e: this->forward_solver->hullMesh->edges()){
    for (BoundaryNormal *tmp_bnd_normal: this->edge_boundary_normals[e]){
      if (tmp_bnd_normal != nullptr && 
          (tmp_bnd_normal->f1 == f || tmp_bnd_normal->f2 == f)){
        Vector3 tmp_normal = tmp_bnd_normal->normal;
        boundary_normals[tmp_bnd_normal->index] = tmp_normal; // + gm_shift;
        // gotta be two neighbors since on edge
        for (BoundaryNormal *neigh_bnd_normal: tmp_bnd_normal->neighbors){ 
          Vector3 neigh_normal = neigh_bnd_normal->normal;
          boundary_normals[neigh_bnd_normal->index] = neigh_normal;
          std::pair<size_t, size_t> tmp_pair = {tmp_bnd_normal->index, neigh_bnd_normal->index};
          if (edge_set.find(tmp_pair) == edge_set.end()){
            edge_set.insert(tmp_pair);
          }
        }
      }
    }
  }
  // set to vector
  std::vector<std::pair<size_t, size_t>> ind_pairs_vector;
  ind_pairs_vector.assign(edge_set.begin(), edge_set.end());
  return {ind_pairs_vector, boundary_normals};
}
        
// ========================================================================================== TOOD: move to optimization and generalize




// hull update stuff
double hull_update_line_search(Eigen::MatrixX3d dfdv, Eigen::MatrixX3d hull_positions, Eigen::Vector3d G_vec, 
                               double bary_reg, double coplanar_reg, double cluster_distance_reg, double unstable_attaction_thresh,
                               std::string policy_general, std::vector<std::pair<Vector3, double>> normal_prob_assignment, 
                               size_t dice_side_count, 
                               double step_size, double decay, bool frozen_G, size_t max_iter, double step_tol){
  
  Forward3DSolver tmp_solver(hull_positions, G_vec, true); // assuming input is convex; will be asserted internally in the constructor
  if (!frozen_G){
    tmp_solver.set_uniform_G();
    G_vec = vec32vec(tmp_solver.get_G());
  }
  tmp_solver.initialize_pre_computes();
  
  // std::cout << "getting fair dice energy for the initial hull\n";
  double min_dice_energy = BoundaryBuilder::dice_energy<double>(hull_positions, G_vec, 
                                                                tmp_solver, bary_reg, coplanar_reg, cluster_distance_reg, unstable_attaction_thresh,
                                                                policy_general, normal_prob_assignment, 
                                                                dice_side_count, false);
  double s = step_size; //

  bool found_smth_optimal = false;
  double tmp_dice_energy;
  int j;
  for (j = 0; j < max_iter; j++) {
    tmp_solver = Forward3DSolver(hull_positions - s * dfdv, G_vec, false); // not necessarily convex
    if (frozen_G && !G_is_inside(*tmp_solver.hullMesh, *tmp_solver.hullGeometry, tmp_solver.get_G())){
        // printf("  - G outside! \n");
        s *= decay;
        continue;
    }
    if (!frozen_G){
        tmp_solver.set_uniform_G();
        G_vec = vec32vec(tmp_solver.get_G());
    }
    tmp_solver.initialize_pre_computes();
    // re-assign since qhull inside solver reshuffles points
    Eigen::MatrixX3d tmp_hull_positions = vertex_data_to_matrix(tmp_solver.hullGeometry->inputVertexPositions);
    Forward3DSolver tmp_solver2(tmp_hull_positions, G_vec, true);
    bool verbose = false;
    if (s < step_tol){
      verbose = true;
      printf("   ---   LS step %d -----\n", j);
    }
    tmp_dice_energy = BoundaryBuilder::dice_energy<double>(tmp_hull_positions, G_vec,
                                                           tmp_solver2, bary_reg, coplanar_reg, cluster_distance_reg, unstable_attaction_thresh,
                                                           policy_general, normal_prob_assignment, 
                                                           dice_side_count, verbose);

    if (tmp_dice_energy < min_dice_energy){
        found_smth_optimal = true;
        break; //  x new is good
    }
    else if(s >= step_tol)
        s *= decay;
    else
      break;
  }
  s = found_smth_optimal ? s : 0.;
  printf("line search for dice ended at iter %d, s: %.10f, \n", j, s);
  // printf("\t\t - line ended at iter %d/%d with s: %f \n", j, max_iter, s);
  return s;
}
