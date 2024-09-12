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

#include "boundary_tools.h"


void draw_arc_on_sphere_static(Vector3 p1, Vector3 p2, Vector3 center, double radius, 
                        size_t seg_count, size_t edge_ind, 
                        double radi_scale, glm::vec3 color, 
                        float arc_curve_radi){
// p1, p2 just represent normal vectors
  if (norm(p1) > 1.01)
    polyscope::warning("p1 norm larger than 1!");
  if (norm(p2) > 1.01)
    polyscope::warning("p2 norm larger than 1!");

  std::vector<std::array<size_t, 2>> edgeInds;
  std::vector<Vector3> positions;
  double sqrt_radi = sqrt(radius);
  // walk on p1-p2 segment
  Vector3 curr_point = p1,
          forward_vec = (p2-p1)/(double)seg_count;
  Vector3 next_point = curr_point + forward_vec;
  Vector3 curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center ,
          next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
  positions.push_back(curr_point_on_sphere);
  for (size_t i = 0; i < seg_count; i++){
    // add to positions list
    curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center ,
    next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
    positions.push_back(next_point_on_sphere);
    // add segment indices
    edgeInds.push_back({i, i+1});

    // update points
    curr_point = next_point;
    next_point += forward_vec;
  }
  auto psArcCurveNet = polyscope::registerCurveNetwork("Arc curve " + std::to_string(edge_ind), positions, edgeInds);
  psArcCurveNet->setRadius(arc_curve_radi * radi_scale, false);
  psArcCurveNet->setColor(color);
  psArcCurveNet->setEnabled(true);
}


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

    // for quick assignment of face-boundary-loops
    face_region_area = FaceData<double>(*forward_solver->hullMesh, 0.);
    // back-flow from all terminal edges
    // printf("  back-flowing terminal edges \n");
    int i = 0;
    for (Edge e: terminal_edges){
        // printf("\n - starting at terminal edge: %d/%d \n", i++, terminal_edges.size());
        // assert(edge_boundary_normals[e].size() == 1); // otherwise we proly have a Gomboc?
        Vector3 stable_edge_normal = forward_solver->edge_stable_normal[e];
        BoundaryNormal *bnd_normal = new BoundaryNormal(stable_edge_normal);
        Face f1 = e.halfedge().face(),
             f2 = e.halfedge().twin().face();
        bnd_normal->f1 = forward_solver->face_last_face[f1];
        bnd_normal->f2 = forward_solver->face_last_face[f2];
        bnd_normal->host_e = e;

        // for visuals
        edge_boundary_normals[e].push_back(bnd_normal);
        
        for (Vertex v: {e.firstVertex(), e.secondVertex()}){
            Vector3 tmp_normal = bnd_normal->normal,
                    f1_normal  = forward_solver->hullGeometry->faceNormal(bnd_normal->f1), // final faces
                    f2_normal  = forward_solver->hullGeometry->faceNormal(bnd_normal->f2),
                    v_normal   = forward_solver->vertex_stable_normal[v];
            // Vector3 imm_f1_normal = forward_solver->hullGeometry->faceNormal(e.halfedge().face()), // immediate face neighbors
            //         imm_f2_normal = forward_solver->hullGeometry->faceNormal(e.halfedge().twin().face());
            
            double f1_area_sign = dot(f1_normal, cross(v_normal, tmp_normal)) >= 0 ? 1. : -1.; // f1 on rhs of bndN->vN
            if (forward_solver->vertex_is_stabilizable[v])
                flow_back_boundary_on_edge(bnd_normal, Edge(), v, 
                                               f1_area_sign);
            else {
                for (Edge neigh_e: v.adjacentEdges()){
                    if (neigh_e != e &&
                        forward_solver->edge_next_vertex[neigh_e] == v){ // neigh_e is a source for this vertex
                        // printf("-SE %d-", neigh_e.getIndex());
                        // face_attraction_boundary[bnd_normal->f1].push_back(bnd_normal);
                        // face_attraction_boundary[bnd_normal->f2].push_back(bnd_normal);
                        flow_back_boundary_on_edge(bnd_normal, neigh_e, v, 
                                                f1_area_sign);
                    }
                }
            }
            // printf(" \n One side done! \n");
        }
    }
}


// recursively follow the boundary curve to a source
bool BoundaryBuilder::flow_back_boundary_on_edge(BoundaryNormal* bnd_normal, Edge src_e, Vertex common_vertex,
                                                 double f1_area_sign){
    // bnd_normal has to be in boundary normals of dest_e; won't assert tho for better performance
    Vertex v = common_vertex; // given as argument for better performance
    // std::cout<<"."<<std::endl;
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
            
            // go with the back-flow
            Vertex next_v = src_e.otherVertex(v);
            if (forward_solver->vertex_is_stabilizable[next_v]){
                // printf(" v ");
                flow_back_boundary_on_edge(new_boundary_normal, Edge(), next_v, f1_area_sign);
            }
            else {
                for (Edge next_src_e: next_v.adjacentEdges()){
                    if (next_src_e != src_e &&
                        forward_solver->edge_next_vertex[next_src_e] == next_v){ // next_src_e is a source for this next vertex
                        // std::cout << " e " << next_src_e.getIndex() << std::endl;
                        // non_singularity and divisive-ness will be checked inside the function
                        bool res = flow_back_boundary_on_edge(new_boundary_normal, next_src_e, next_v, f1_area_sign);
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
    double curr_f1_alignment = dot(f1_normal, cross(next_normal, tmp_normal)) >= 0 ? 1. : -1.; // checking alignment again since it could change along the way
    double curr_f2_alignment = dot(f2_normal, cross(next_normal, tmp_normal)) >= 0 ? 1. : -1.;
    double f1_sign_change = f1_area_sign == curr_f1_alignment ? 1. : -1;
    double f2_sign_change = (-f1_area_sign == curr_f2_alignment) ? 1.: -1;
    face_region_area[bnd_normal->f1] += f1_sign_change * abs(triangle_patch_area_on_sphere(f1_normal, bnd_normal->normal, next_normal));
    face_region_area[bnd_normal->f2] += f2_sign_change * abs(triangle_patch_area_on_sphere(f2_normal, bnd_normal->normal, next_normal));
    // printf(" -ret- ");
    return true;
    // TODO: take care of when f1,f2 on the same side when starting from saddle 
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
    for (Face f: forward_solver->hullMesh->faces()){
        if (face_region_area[f] > 0){
            // face_region_area[f] /= (4.*PI);
            probs.push_back({f, face_region_area[f]/(4.*PI)});
        }
    }
    std::sort(probs.begin(), probs.end(), [] (auto a, auto b) { return a.second > b.second; });
    printf("sorted probs: \n");
    for (auto pair: probs)
        printf("  -f %d: %f\n", pair.first.getIndex(), pair.second);
}









// autodiff stuff
// --------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------









void BoundaryBuilder::build_boundary_normals_for_autodiff( // autodiff::MatrixX3var &var_positions, autodiff::Vector3var &var_G,
                                                          bool generate_gradients){
    
    // assigning the autodiff variables
    var_positions = vertex_data_to_matrix(forward_solver->hullGeometry->inputVertexPositions);
    var_positions_vec = autodiff::VectorXvar{var_positions.reshaped()};
    // var_positions = VertexData<autodiff::Vector3var>(*forward_solver->hullMesh);
    // for (Vertex v: forward_solver->hullMesh->vertices())
    //     var_positions.row(v.getIndex()) = vec32vec(forward_solver->hullGeometry->inputVertexPositions[v]);
    var_G = autodiff::Vector3var{forward_solver->get_G().x, forward_solver->get_G().y, forward_solver->get_G().z};//vec32vec(forward_solver->get_G());

    //
    face_normals_ad = FaceData<autodiff::Vector3var>(*forward_solver->hullMesh);
    for (Face f: forward_solver->hullMesh->faces()){
        Halfedge he = f.halfedge();
        Vertex v1 = he.vertex(),
               v2 = he.next().vertex(),
               v3 = he.next().next().vertex();
        // face_normals_ad[f] = (var_positions.row(v2.getIndex()) - var_positions.row(v1.getIndex())).cross(var_positions.row(v3.getIndex()) - var_positions.row(v1.getIndex())).normalized();
        face_normals_ad[f] = (var_positions.row(v2.getIndex()) - var_positions.row(v1.getIndex())).cross(var_positions.row(v3.getIndex()) - var_positions.row(v1.getIndex()));//.normalized();
    }

    //
    vertex_boundary_normal = VertexData<BoundaryNormal*>(*forward_solver->hullMesh, nullptr);
    edge_boundary_normals  = EdgeData<std::vector<BoundaryNormal*>>(*forward_solver->hullMesh);
    // edge_boundary_normals = EdgeData<std::vector<Vector3>>(*forward_solver->hullMesh);

    if (generate_gradients){
        Eigen::MatrixX3d zero_mat = Eigen::MatrixX3d::Zero(var_positions.rows(), 3);
        df_dv_grads_ad = FaceData<Eigen::MatrixX3d>(*forward_solver->hullMesh, zero_mat);
        Eigen::Vector3d zero_vec = Eigen::Vector3d::Zero();
        df_dG_grads = FaceData<Eigen::Vector3d>(*forward_solver->hullMesh, zero_vec);
    }
    BoundaryNormal::counter = 0;
    // 
    printf("  finding terminal edges \n");

    std::vector<Edge> terminal_edges;
    for (Edge e: forward_solver->hullMesh->edges()){
        if (forward_solver->edge_next_vertex[e].getIndex() == INVALID_IND){ // singular edge
            Face f1 = e.halfedge().face(),
                 f2 = e.halfedge().twin().face();
            if (forward_solver->face_last_face[f1] != forward_solver->face_last_face[f2]){ // saddle edge
                // proved: at least one singular edge like this must exist
                // TODO: assert that stable normal falls inside the edge arc 
                terminal_edges.push_back(e);
            }
        }
    }

    // for quick assignment of face-boundary-loops
    face_region_area = FaceData<double>(*forward_solver->hullMesh, 0.);
    face_region_area_complex_value_ad = FaceData<autodiff::Vector2var>(*forward_solver->hullMesh, autodiff::Vector2var(1., 0.));
    // face_region_area_ad = FaceData<autodiff::var>(*forward_solver->hullMesh, 0.);
    // back-flow from all terminal edges
    printf("  back-flowing terminal edges %d \n", terminal_edges.size());
    int i = 0;
    for (Edge e: terminal_edges){
        printf("\n - starting at terminal edge: %d/%d \n", i++, terminal_edges.size());
        // assert(edge_boundary_normals[e].size() == 1); // otherwise we proly have a Gomboc!
        Vector3 stable_edge_normal = forward_solver->edge_stable_normal[e];
        BoundaryNormal *bnd_normal = new BoundaryNormal(stable_edge_normal);
        // bnd_normal->normal_ad = point_to_segment_normal_ad(e);
        autodiff::Vector3var bnd_normal_ad = point_to_segment_normal_ad(e);
       
        Face f1 = e.halfedge().face(),
             f2 = e.halfedge().twin().face();
        bnd_normal->f1 = forward_solver->face_last_face[f1];
        bnd_normal->f2 = forward_solver->face_last_face[f2];
        bnd_normal->host_e = e;
        // for visuals
        edge_boundary_normals[e].push_back(bnd_normal);
        for (Vertex v: {e.firstVertex(), e.secondVertex()}){
            
            // // TESTING dependency speed-ups
            // std::set<Vertex> effective_vertices;
            // effective_vertices.insert(e.firstVertex());
            // effective_vertices.insert(e.secondVertex());
            // for (Vertex tmp_v: bnd_normal->f1.adjacentVertices())
            //     effective_vertices.insert(tmp_v);
            // for (Vertex tmp_v: bnd_normal->f2.adjacentVertices())
            //     effective_vertices.insert(tmp_v);
            
            Vector3 tmp_normal = bnd_normal->normal,
                    f1_normal  = forward_solver->hullGeometry->faceNormal(bnd_normal->f1),
                    f2_normal  = forward_solver->hullGeometry->faceNormal(bnd_normal->f2),
                    v_normal   = forward_solver->vertex_stable_normal[v];
            // Vector3 imm_f1_normal = forward_solver->hullGeometry->faceNormal(e.halfedge().face()), // immediate face neighbors
            //         imm_f2_normal = forward_solver->hullGeometry->faceNormal(e.halfedge().twin().face());
            double f1_area_sign = dot(f1_normal, cross(v_normal, tmp_normal)) >= 0 ? 1. : -1.; // f1 on rhs of bndN->vN
    
            if (forward_solver->vertex_is_stabilizable[v]){
                // printf("-SV-\n");
                flow_back_boundary_on_edge_for_autodiff(bnd_normal, bnd_normal_ad, Edge(), v, f1_area_sign //, var_positions, var_G
                                                        // ,effective_vertices
                                                        );
            }
            else {
                for (Edge neigh_e: v.adjacentEdges()){
                    if (neigh_e != e &&
                        forward_solver->edge_next_vertex[neigh_e] == v){ // neigh_e is a source for this vertex
                        // face_attraction_boundary[bnd_normal->f1].push_back(bnd_normal);
                        // face_attraction_boundary[bnd_normal->f2].push_back(bnd_normal);
                        // printf(" -Se %d-", neigh_e.getIndex());
                        bool res = flow_back_boundary_on_edge_for_autodiff(bnd_normal, bnd_normal_ad, neigh_e, v, f1_area_sign //, effective_vertices //, var_positions, var_G
                        );
                        if (res) break; // kinda redundant; since source-ness is checked before going in
                    }
                }
            }
            printf(" \n One side done! \n");
        }
    }
    
    // double total_area = 0.;
    // autodiff::var total_area_ad = 0.;
    // for (Face f: forward_solver->hullMesh->faces()){
    //     if (forward_solver->face_last_face[f] == f){
    //         total_area += face_region_area[f];
    //         total_area_ad += 2 * atan2(face_region_area_complex_value_ad[f][1], face_region_area_complex_value_ad[f][0]);
    //         // std::cout << "face " << f.getIndex() << " area: " << face_region_area[f] << ",  " << face_region_area_ad[f] << std::endl;
    //     }
    // }
    // std::cout << "total face areas: " << total_area << "--- from complex stuff: " << total_area_ad << std::endl;
    
    // -------- gradients --------
    // if (generate_gradients){
    //     // df/dv
    //     for (Face f: forward_solver->hullMesh->faces()){
    //         if (forward_solver->face_last_face[f] == f){
    //             autodiff::Vector2var polygon_complex_ad = face_region_area_complex_value_ad[f];
    //             autodiff::var polygon_area = 2 * atan2(polygon_complex_ad[1], polygon_complex_ad[0]);
    //             // autodiff::VectorXvar poses_ad_vec = autodiff::VectorXvar{var_positions.reshaped()};
    //             Eigen::VectorXd dfdv = autodiff::gradient(polygon_area, var_positions_vec);
    //             Eigen::MatrixXd dfdv_mat = dfdv.reshaped(var_positions.rows(), var_positions.cols());
    //             df_dv_grads_ad[f] = dfdv_mat;
    //             // df/dG
    //             df_dG_grads[f] = autodiff::gradient(polygon_area, var_G).cast<double>();
    //         }
    //     }
    // }
}


bool BoundaryBuilder::flow_back_boundary_on_edge_for_autodiff(BoundaryNormal* bnd_normal, autodiff::Vector3var bnd_normal_ad, Edge src_e, 
                                                        Vertex common_vertex, double f1_area_sign
                                                        //, std::set<Vertex> &effective_vertices
                                                        //, autodiff::MatrixX3var &var_positions, autodiff::Vector3var &var_G
                                                        ){
    // bnd_normal has to be in boundary normals of dest_e; won't assert tho for better performance
    Vertex v = common_vertex; // given as argument for better performance
    // std::cout<<" . "<<std::endl;
    
    BoundaryNormal* next_bnd_normal = nullptr;
    autodiff::Vector3var next_bnd_normal_ad;
    if (forward_solver->vertex_is_stabilizable[v]){ // we are at a source
        // BoundaryNormal *curr_vertex_boundary_normal = vertex_boundary_normal[v];
        next_bnd_normal = vertex_boundary_normal[v];
        if (vertex_boundary_normal[v] == nullptr){ // first time arriving at this stable vertex; create the boundary normal
            Vector3 stable_vertex_normal = forward_solver->vertex_stable_normal[v];
            next_bnd_normal = new BoundaryNormal(stable_vertex_normal);
            // ad stuff
            // next_bnd_normal->normal_ad = var_positions.row(v.getIndex()) - var_G.transpose(); // not normalizing; for more stability of gradients?
            
            vertex_boundary_normal[v] = next_bnd_normal;
            next_bnd_normal->host_v = v;
        }
        next_bnd_normal_ad = (var_positions.row(v.getIndex()) - var_G.transpose());//.normalized(); // var_G
        next_bnd_normal->add_neighbor(bnd_normal);
        bnd_normal->add_neighbor(next_bnd_normal);
        // effective_vertices.insert(v);
        evaluate_boundary_patch_area_and_grads_ad(bnd_normal, bnd_normal_ad, 
                                                  next_bnd_normal, next_bnd_normal_ad, 
                                                  f1_area_sign  //   ,effective_vertices
                                                  );
        printf("  -> V done! \n");
    }
    else {
        // vertex is not an equilibria; i.e. source is outside the vertex
        if (forward_solver->edge_next_vertex[src_e] == v){ // src_e is a source for this vertex
            Face f1 = src_e.halfedge().face(),
                 f2 = src_e.halfedge().twin().face();
            Vector3 vertex_stable_normal = forward_solver->vertex_stable_normal[v];
            bool sign_change = false;
            Vector3 e_bnd_normal = intersect_arc_ray_with_arc(vertex_stable_normal, bnd_normal->normal, 
                                                              forward_solver->hullGeometry->faceNormal(f1), 
                                                              forward_solver->hullGeometry->faceNormal(f2), sign_change);
            if(e_bnd_normal.norm() == 0.){
                // printf(" -src but miss- ");
                return false; // not a source for the given bnd_normal
            }
            printf(" -> e %d ", src_e.getIndex());
            // found the source normal
            // BoundaryNormal *new_boundary_normal = new BoundaryNormal(e_bnd_normal);
            next_bnd_normal = new BoundaryNormal(e_bnd_normal);
            edge_boundary_normals[src_e].push_back(next_bnd_normal);

            // ad stuff
            next_bnd_normal_ad = intersect_arcs_ad(v, bnd_normal_ad, face_normals_ad[f1], face_normals_ad[f2], sign_change);
            // next_bnd_normal->normal_ad = intersect_arcs_ad(v, bnd_normal->normal_ad, face_normals_ad[f1], face_normals_ad[f2], sign_change);
            
            // intersect_arc_ray_with_arc_ad(var_positions, var_G, v, bnd_normal->normal_ad, tmp_f1_normal_ad, tmp_f2_normal_ad, sign_change);

            bnd_normal->add_neighbor(next_bnd_normal);
            next_bnd_normal->add_neighbor(bnd_normal);
            next_bnd_normal->host_e = src_e;

            // region labels
            next_bnd_normal->f1 = bnd_normal->f1;
            next_bnd_normal->f2 = bnd_normal->f2;
            
            printf(" eval \n");
            // effective_vertices.insert(v);
            // for (Vertex tmp_v: f1.adjacentVertices())
            //     effective_vertices.insert(tmp_v);
            // for (Vertex tmp_v: f2.adjacentVertices())
            //     effective_vertices.insert(tmp_v);

            evaluate_boundary_patch_area_and_grads_ad(bnd_normal, bnd_normal_ad, 
                                                      next_bnd_normal, next_bnd_normal_ad, 
                                                      f1_area_sign
                                                    //   ,effective_vertices
                                                      );
            // go with the back-flow
            // printf("recursion back flow!\n");
            Vertex next_v = src_e.otherVertex(v);
            if (forward_solver->vertex_is_stabilizable[next_v]){
                printf("v ");
                flow_back_boundary_on_edge_for_autodiff(next_bnd_normal, next_bnd_normal_ad, Edge(), next_v, f1_area_sign //, var_positions, var_G
                                                        // ,effective_vertices
                                                        );
            }
            else{
                for (Edge next_src_e: next_v.adjacentEdges()){
                    if (next_src_e != src_e && 
                        forward_solver->edge_next_vertex[next_src_e] == next_v){
                        // non_singularity and divisive-ness will be checked inside the function
                        // std::cout << " e " << next_src_e.getIndex();
                        bool res = flow_back_boundary_on_edge_for_autodiff(next_bnd_normal, next_bnd_normal_ad, next_src_e, next_v, f1_area_sign //, var_positions, var_G
                                                                        //    ,effective_vertices
                                                                           );
                        if (res) // got to maximum and returning
                            break;
                    }
                }
            }
        }
    }
    if (bnd_normal->normal.norm() == 0. || next_bnd_normal == nullptr) { // not src
        // std::cout << "ret \n" << std::endl;
        return false;
    }
    printf(" -ret- \n");
    return true;
    // TODO: take care of when f1,f2 on the same side when starting from saddle 
}


void BoundaryBuilder::evaluate_boundary_patch_area_and_grads_ad(BoundaryNormal* bnd_normal1, autodiff::Vector3var bnd_normal1_ad, 
                                                                BoundaryNormal* bnd_normal2, autodiff::Vector3var bnd_normal2_ad, 
                                                                double f1_area_sign,
                                                                bool complex_accum
                                                                // , std::set<Vertex> &effective_vertices
                                                                ){
    // printf("None-ad normals and patches..\n");
    Vector3 f1_normal = forward_solver->hullGeometry->faceNormal(bnd_normal1->f1), // f1,f2 are the same along the current path to maximum
            f2_normal = forward_solver->hullGeometry->faceNormal(bnd_normal1->f2);
    // TODO: flipped after static function correction
    double curr_f1_alignment = dot(f1_normal, cross(bnd_normal1->normal, bnd_normal2->normal)) >= 0 ? 1. : -1.; // checking alignment again since it could change along the way
    double curr_f2_alignment = dot(f2_normal, cross(bnd_normal1->normal, bnd_normal2->normal)) >= 0 ? 1. : -1.;
    double f1_sign_change = f1_area_sign == curr_f1_alignment ? 1. : -1;
    double f2_sign_change = (-f1_area_sign == curr_f2_alignment) ? 1.: -1;
    double f1_side_patch = triangle_patch_area_on_sphere(f1_normal, bnd_normal1->normal, bnd_normal2->normal),
           f2_side_patch = triangle_patch_area_on_sphere(f2_normal, bnd_normal1->normal, bnd_normal2->normal);
    bool flip_f1 = f1_side_patch < 0.,
         flip_f2 = f2_side_patch < 0.;
    face_region_area[bnd_normal1->f1] += f1_sign_change * f1_side_patch; // TODO : removing abs with new sign
    face_region_area[bnd_normal1->f2] += f2_sign_change * f2_side_patch;
    // printf(" getting AD patches ..\n");
    // DEBUG
    
    // TESTING speed up

    if (complex_accum){
        autodiff::Vector2var f1_side_patch_ad_complex = solid_angle_complex_ad(face_normals_ad[bnd_normal1->f1], bnd_normal1_ad, bnd_normal2_ad),
                             f2_side_patch_ad_complex = solid_angle_complex_ad(face_normals_ad[bnd_normal1->f2], bnd_normal1_ad, bnd_normal2_ad);
        if (flip_f1)
            f1_side_patch_ad_complex[1] = -1. * f1_side_patch_ad_complex[1];
        if (f1_sign_change == -1)
            f1_side_patch_ad_complex[1] = -1. * f1_side_patch_ad_complex[1];
        if (flip_f2)
            f2_side_patch_ad_complex[1] = -1. * f2_side_patch_ad_complex[1];
        if (f2_sign_change == -1)
            f2_side_patch_ad_complex[1] = -1. * f2_side_patch_ad_complex[1];
        face_region_area_complex_value_ad[bnd_normal1->f1] = complex_mult(face_region_area_complex_value_ad[bnd_normal1->f1],
                                                                        f1_side_patch_ad_complex);
        face_region_area_complex_value_ad[bnd_normal1->f2] = complex_mult(face_region_area_complex_value_ad[bnd_normal1->f2],
                                                                        f2_side_patch_ad_complex);
    }
    else { // Without complex stuff
        autodiff::var f1_side_patch_ad = triangle_patch_area_on_sphere_ad(face_normals_ad[bnd_normal1->f1], bnd_normal1_ad, bnd_normal2_ad),
                      f2_side_patch_ad = triangle_patch_area_on_sphere_ad(face_normals_ad[bnd_normal1->f2], bnd_normal1_ad, bnd_normal2_ad);
        // f1
        // printf("grads ..\n");
        if (f1_side_patch != 0.){
            Eigen::VectorXd df1dv = autodiff::gradient(f1_side_patch_ad, var_positions_vec);
            Eigen::MatrixXd df1dv_mat = df1dv.reshaped(var_positions.rows(), var_positions.cols());
            df_dv_grads_ad[bnd_normal1->f1] += curr_f1_alignment * df1dv_mat;
            // printf(" dfdv norm: %f\n", df_dv_grads_ad[bnd_normal1->f1].norm());
            // for (Vertex v: effective_vertices){
            //     df_dv_grads_ad[bnd_normal1->f1][v] += autodiff::gradient(f1_side_patch_ad, var_positions.row(v.getIndex()));
            // }
            // df/dG
            df_dG_grads[bnd_normal1->f1] += curr_f1_alignment * autodiff::gradient(f1_side_patch_ad, var_G);
        }
        // f2
        if (f2_side_patch != 0.){
            Eigen::VectorXd df2dv = autodiff::gradient(f2_side_patch_ad, var_positions_vec);
            Eigen::MatrixXd df2dv_mat = df2dv.reshaped(var_positions.rows(), var_positions.cols());
            df_dv_grads_ad[bnd_normal1->f2] += curr_f2_alignment * df2dv_mat;
            // for (Vertex v: effective_vertices){
            //     df_dv_grads_ad[bnd_normal1->f2][v] += autodiff::gradient(f2_side_patch_ad, var_positions.row(v.getIndex()));
            // }
            // df/dG
            df_dG_grads[bnd_normal1->f2] += curr_f2_alignment * autodiff::gradient(f2_side_patch_ad, var_G);    
        }
    }
    // }
    // if (f1_side_patch != 0.) // avoiding conditions in ad version since idk how they work yet
    //     face_region_area_ad[bnd_normal->f1] += f1_sign_change * f1_side_patch_ad;
    // if (f2_side_patch != 0.)
    //     face_region_area_ad[bnd_normal->f2] += f2_sign_change * f2_side_patch_ad;
    // if (verbose){
    //     printf("  - NO AD %f, %f\n", face_region_area[bnd_normal->f1], face_region_area[bnd_normal->f2]);
    //     std::cout << "  -    AD " << face_region_area_ad[bnd_normal->f1] << ", " << face_region_area_ad[bnd_normal->f2] << "\n";
    //     std::cout << "  - patches   " << f1_side_patch << ", " << f2_side_patch << "\n";
    //     std::cout << "  - AD patchs " << f1_side_patch_ad << ", " << f2_side_patch_ad << "\n";
    //     printf("--------------------------------------------\n");
    // }
    // printf("area compute done!\n");
}

autodiff::Vector3var BoundaryBuilder::point_to_segment_normal_ad(Edge e){
    Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
    autodiff::Vector3var PB = var_positions.row(v2.getIndex()) - var_G.transpose(),// var_G,
                         AB = var_positions.row(v2.getIndex()) - var_positions.row(v1.getIndex());
    return (PB * AB.dot(AB) - AB * AB.dot(PB));//.normalized();
}

autodiff::Vector3var BoundaryBuilder::intersect_arcs_ad(Vertex v, autodiff::Vector3var &R2, autodiff::Vector3var &A, autodiff::Vector3var &B, bool sign_change){
    autodiff::Vector3var p = (((var_positions.row(v.getIndex()) - var_G.transpose()).cross(R2)).cross(A.cross(B)));//.normalized(); // var_G.transpose()
    if (sign_change)
        return -p;
    return p;
}
        
autodiff::var BoundaryBuilder::triangle_patch_area_on_sphere_ad(autodiff::Vector3var &A, autodiff::Vector3var &B, autodiff::Vector3var &C){
    // autodiff::Vector3var n_AB = A.cross(B),
    //                      n_BC = B.cross(C),
    //                      n_CA = C.cross(A);
    // autodiff::var angle_A = atan2(n_AB.cross(-n_CA).norm(), n_AB.dot(-n_CA)),
                //   angle_B = atan2(n_BC.cross(-n_AB).norm(), n_BC.dot(-n_AB)),
    //               angle_C = atan2(n_CA.cross(-n_BC).norm(), n_CA.dot(-n_BC));
    // return angle_A + angle_B + angle_C - PI;
    autodiff::var Adotbc = A.dot(B.cross(C));
    autodiff::var denom = A.norm()*B.norm()*C.norm() + A.dot(B)*C.norm() + B.dot(C)*A.norm() + C.dot(A)*B.norm();
    return 2. * atan2(Adotbc, denom);
}


autodiff::Vector2var BoundaryBuilder::solid_angle_complex_ad(autodiff::Vector3var &A, autodiff::Vector3var &B, autodiff::Vector3var &C){
    autodiff::var Adotbc = A.dot(B.cross(C));
    autodiff::var denom = A.norm()*B.norm()*C.norm() + A.dot(B)*C.norm() + B.dot(C)*A.norm() + C.dot(A)*B.norm();
    return autodiff::Vector2var{denom, Adotbc};
}

autodiff::Vector2var BoundaryBuilder::complex_mult(autodiff::Vector2var &A, autodiff::Vector2var &B){
    return autodiff::Vector2var{A[0]*B[0] - A[1]*B[1], A[0]*B[1] + A[1]*B[0]};
}



// ========================================================================================== TOOD: move to optimization and generalize




// hull update stuff
// frozen G so far
double hull_update_line_search(Eigen::MatrixX3d dfdv, Eigen::MatrixX3d hull_positions, Eigen::Vector3d G_vec,
                               size_t dice_side_count, 
                               double step_size, double decay, bool frozen_G, size_t max_iter){
  
  Forward3DSolver tmp_solver(hull_positions, G_vec);
  tmp_solver.initialize_pre_computes();
  BoundaryBuilder tmp_builder(&tmp_solver);
  tmp_builder.build_boundary_normals();

  double min_dice_energy = BoundaryBuilder::dice_energy<double>(hull_positions, G_vec, tmp_solver, dice_side_count);
//   printf(" current fair dice energy: %f\n", s_min_dice_energy);
  double s = step_size; //

  bool found_smth_optimal = false;
  int j;
  double tmp_dice_energy;
  for (j = 0; j < max_iter; j++) {
      // update stuff
    //   printf(" ^^ at line search iter: %d  s = %f, DE: %f\n", j, s, tmp_fair_dice_energy);
    //   auto tmpmesh = polyscope::registerSurfaceMesh("temp hull LINE SEARCH", initial_hull_poses + s * grad,
    //                                                 fwd_solver.hullMesh->getFaceVertexList());
    //   tmpmesh->setSurfaceColor({0.6,0,0.5});
    //   polyscope::show();
        tmp_solver = Forward3DSolver(hull_positions + s * dfdv, G_vec);
        if (!G_is_inside(*tmp_solver.hullMesh, *tmp_solver.hullGeometry, tmp_solver.get_G())){
            printf("  - G outside! \n");
            s *= decay;
            continue;
        }
        tmp_solver.initialize_pre_computes();
        tmp_dice_energy = BoundaryBuilder::dice_energy<double>(hull_positions + s * dfdv,
                             G_vec, tmp_solver, dice_side_count);
    //   printf("  *** temp fair dice energy %d: %f\n", j, tmp_fair_dice_energy);

      if (tmp_dice_energy < min_dice_energy){
        found_smth_optimal = true;
        break; //  x new is good
      }
      else
          s *= decay;
  }
  s = found_smth_optimal ? s : 0.;
//   printf("line search for dice ended at iter %d, s: %.10f, \n \t\t\t\t\t fnew: %f \n", j, s, tmp_fair_dice_energy);
  printf("\t\t - line ended at iter %d/%d\n", j, max_iter);
  return s;
}

// ========================================================================================== STATIC FUNCTIONS



