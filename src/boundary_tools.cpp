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


size_t BoundaryNormal::counter = 0;

// constructor
BoundaryNormal::BoundaryNormal(Vector3 _normal)
    : index(counter++){
        normal = _normal;
        host_e = Edge();
        host_v = Vertex();
}

// void BoundaryNormal::add_neighbor(BoundaryNormal* _neigh) {
//     neighbors.push_back(_neigh);
// }

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
    return terminal_edges;
}

void BoundaryBuilder::build_boundary_normals(){
    vertex_boundary_normal = VertexData<BoundaryNormal*>(*forward_solver->hullMesh, nullptr);
    // edge_boundary_normals  = EdgeData<std::vector<BoundaryNormal*>>(*forward_solver->hullMesh);
    edge_boundary_normals = EdgeData<std::vector<Vector3>>(*forward_solver->hullMesh);
    BoundaryNormal::counter = 0;
    // 
    // printf("  buidling face-last-face\n");
    forward_solver->build_face_last_faces(); // calls face next face within
    // printf("  finding terminal edges \n");
    std::vector<Edge> terminal_edges = find_terminal_edges();

    // for quick assignment of face-boundary-loops
    face_region_area = FaceData<double>(*forward_solver->hullMesh, 0.);
    // back-flow from all terminal edges
    // printf("  back-flowing terminal edges \n");
    int i = 0;
    for (Edge e: terminal_edges){
        printf("\n - starting at terminal edge: %d/%d \n", i++, terminal_edges.size());
        // assert(edge_boundary_normals[e].size() == 1); // otherwise we proly have a Gomboc?
        Vector3 stable_edge_normal = forward_solver->edge_stable_normal[e];
        BoundaryNormal *bnd_normal = new BoundaryNormal(stable_edge_normal);
        Face f1 = e.halfedge().face(),
             f2 = e.halfedge().twin().face();
        bnd_normal->f1 = forward_solver->face_last_face[f1];
        bnd_normal->f2 = forward_solver->face_last_face[f2];
        bnd_normal->host_e = e;

        // for visuals
        edge_boundary_normals[e].push_back(bnd_normal->normal);
        
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
                    if (neigh_e != e){
                        printf("-SE %d-", neigh_e.getIndex());
                        // face_attraction_boundary[bnd_normal->f1].push_back(bnd_normal);
                        // face_attraction_boundary[bnd_normal->f2].push_back(bnd_normal);
                        flow_back_boundary_on_edge(bnd_normal, neigh_e, v, 
                                                f1_area_sign);
                    }
                }
            }
        }
    }
}




// recursively follow the boundary curve to a source
bool BoundaryBuilder::flow_back_boundary_on_edge(BoundaryNormal* bnd_normal, Edge src_e, Vertex common_vertex,
                                                 double f1_area_sign){
    // bnd_normal has to be in boundary normals of dest_e; won't assert tho for better performance
    Vertex v = common_vertex; // given as argument for better performance
    std::cout<<"."<<std::flush;
    Vector3 tmp_normal({0.,0.,0.}), next_normal({0.,0.,0.});
    if (forward_solver->vertex_is_stabilizable[v]){ // we are at a source
        BoundaryNormal *curr_vertex_boundary_normal = vertex_boundary_normal[v];
        if (vertex_boundary_normal[v] == nullptr){ // first time arriving at this stable vertex; create the boundary normal
            Vector3 stable_vertex_normal = forward_solver->vertex_stable_normal[v];
            curr_vertex_boundary_normal = new BoundaryNormal(stable_vertex_normal);
            vertex_boundary_normal[v] = curr_vertex_boundary_normal;
            curr_vertex_boundary_normal->host_v = v;
        }
        // curr_vertex_boundary_normal->add_neighbor(bnd_normal);
        // bnd_normal->add_neighbor(curr_vertex_boundary_normal);
        // handling face region boundary
        // face_attraction_boundary[bnd_normal->f1].push_back(curr_vertex_boundary_normal);
        // face_attraction_boundary[bnd_normal->f2].push_back(curr_vertex_boundary_normal);

        tmp_normal = bnd_normal->normal;
        next_normal = curr_vertex_boundary_normal->normal;
        std::cout<<" V! "<<std::flush;

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
                printf(" -src but miss- ");
                return false; // not a source for the given bnd_normal
            }
            // found the source normal


            BoundaryNormal *new_boundary_normal = new BoundaryNormal(e_bnd_normal);
            // edge_boundary_normals[src_e].push_back(new_boundary_normal);
            edge_boundary_normals[src_e].push_back(e_bnd_normal);

            // bnd_normal->add_neighbor(new_boundary_normal);
            // new_boundary_normal->add_neighbor(bnd_normal);
            new_boundary_normal->host_e = src_e;

            // region labels
            new_boundary_normal->f1 = bnd_normal->f1;
            new_boundary_normal->f2 = bnd_normal->f2;

            tmp_normal  = bnd_normal->normal,
            next_normal = new_boundary_normal->normal;
            
            // go with the back-flow
            Vertex next_v = src_e.otherVertex(v);
            if (forward_solver->vertex_is_stabilizable[next_v]){
                std::cout << "v" << std::flush;           
                flow_back_boundary_on_edge(new_boundary_normal, Edge(), next_v, f1_area_sign);
            }
            else{
                for (Edge next_src_e: next_v.adjacentEdges()){
                    if (next_src_e != src_e){
                        std::cout << "e " << next_src_e.getIndex() << std::flush;
                        // non_singularity and divisive-ness will be checked inside the function
                        bool res = flow_back_boundary_on_edge(new_boundary_normal, next_src_e, next_v, f1_area_sign);
                        if (res) // got to maximum and returning
                            break;
                    }
                }
            }
        }
        else
            std::cout << " not src! " << std::flush;
    }
    std::cout << "F" << std::flush;
    
    if (tmp_normal.norm() == 0.){ // when does this happen?
        std::cout << "ret\n" << std::flush;
        return false;
    }
    Vector3 f1_normal = forward_solver->hullGeometry->faceNormal(bnd_normal->f1), // f1,f2 are the same along the current path to maximum
            f2_normal = forward_solver->hullGeometry->faceNormal(bnd_normal->f2);
    double curr_f1_alignment = dot(f1_normal, cross(next_normal, tmp_normal)) >= 0 ? 1. : -1.; // checking alignment again since it could change along the way
    double curr_f2_alignment = dot(f2_normal, cross(next_normal, tmp_normal)) >= 0 ? 1. : -1.;
    double f1_sign_change = f1_area_sign == curr_f1_alignment ? 1. : -1;
    double f2_sign_change = (-f1_area_sign == curr_f2_alignment) ? 1.: -1;
    if (f1_sign_change == 1)
        face_region_area[bnd_normal->f1] += triangle_patch_area_on_sphere(f1_normal, bnd_normal->normal, next_normal);
    else
        face_region_area[bnd_normal->f1] -= triangle_patch_area_on_sphere(f1_normal, bnd_normal->normal, next_normal);
    if (f2_sign_change == 1)
        face_region_area[bnd_normal->f2] += triangle_patch_area_on_sphere(f2_normal, bnd_normal->normal, next_normal);
    else 
        face_region_area[bnd_normal->f2] -= triangle_patch_area_on_sphere(f2_normal, bnd_normal->normal, next_normal);
    
    printf(" gg \n");
    return true;
    // TODO: take care of when f1,f2 on the same side when starting from saddle 
}


double BoundaryBuilder::get_fair_dice_energy(size_t side_count){
    double energy = 0., goal_area = 4.*PI/(double)side_count;
    std::vector<double> face_areas;
    for (Face f: forward_solver->hullMesh->faces()){
        if (forward_solver->face_last_face[f] == f)
            face_areas.push_back(face_region_area[f]);
    }
    std::sort(face_areas.begin(), face_areas.end());
    double nth_largest_area = face_areas[face_areas.size() - std::min(side_count, face_areas.size())];
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
    return energy;
}

void BoundaryBuilder::print_area_of_boundary_loops(){
    printf(" Face probs:\n");
    std::vector<double> probs;
    for (Face f: forward_solver->hullMesh->faces()){
        if (face_region_area[f] > 0){
            // face_region_area[f] /= (4.*PI);
            probs.push_back(face_region_area[f]/(4.*PI));
            // printf(" f %d: %f\n", f.getIndex(), face_region_area[f]/(4.*PI));
        }
    }
    std::sort(probs.begin(), probs.end(), [] (auto a, auto b) { return a > b; });
    printf("sorted probs: \n");
    for (double prob: probs)
        printf("  -%f\n", prob);
}




// autodiff stuff
// --------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------




void BoundaryBuilder::build_boundary_normals_for_autodiff(autodiff::MatrixX3var &var_positions, autodiff::Vector3var &var_G,
                                                          bool generate_gradients){
    vertex_boundary_normal = VertexData<BoundaryNormal*>(*forward_solver->hullMesh, nullptr);
    // edge_boundary_normals  = EdgeData<std::vector<BoundaryNormal*>>(*forward_solver->hullMesh);
    edge_boundary_normals = EdgeData<std::vector<Vector3>>(*forward_solver->hullMesh);

    if (generate_gradients){
        Eigen::MatrixX3d zero_mat = Eigen::MatrixX3d::Zero(var_positions.rows(), 3);
        df_dv_grads_ad = FaceData<Eigen::MatrixX3d>(*forward_solver->hullMesh, zero_mat);
        Eigen::Vector3d zero_vec = Eigen::Vector3d::Zero();
        df_dG_grads = FaceData<Eigen::Vector3d>(*forward_solver->hullMesh, zero_vec);
    }
    BoundaryNormal::counter = 0;
    // 
    // printf("  buidling face-last-face\n");
    forward_solver->build_face_last_faces(); // calls face next face within
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
    face_region_area_ad = FaceData<autodiff::var>(*forward_solver->hullMesh, 0.);
    // back-flow from all terminal edges
    printf("  back-flowing terminal edges \n");
    int i = 0;
    for (Edge e: terminal_edges){
        printf("\n - starting at terminal edge: %d/%d \n", i++, terminal_edges.size());
        // assert(edge_boundary_normals[e].size() == 1); // otherwise we proly have a Gomboc!
        Vector3 stable_edge_normal = forward_solver->edge_stable_normal[e];
        BoundaryNormal *bnd_normal = new BoundaryNormal(stable_edge_normal);
        bnd_normal->normal_ad = point_to_segment_normal_ad(var_positions, var_G, e);
        Face f1 = e.halfedge().face(),
             f2 = e.halfedge().twin().face();
        bnd_normal->f1 = forward_solver->face_last_face[f1];
        bnd_normal->f2 = forward_solver->face_last_face[f2];
        bnd_normal->host_e = e;
        // for visuals
        edge_boundary_normals[e].push_back(bnd_normal->normal);
        for (Vertex v: {e.firstVertex(), e.secondVertex()}){
            Vector3 tmp_normal = bnd_normal->normal,
                    f1_normal  = forward_solver->hullGeometry->faceNormal(bnd_normal->f1),
                    f2_normal  = forward_solver->hullGeometry->faceNormal(bnd_normal->f2),
                    v_normal   = forward_solver->vertex_stable_normal[v];
            // Vector3 imm_f1_normal = forward_solver->hullGeometry->faceNormal(e.halfedge().face()), // immediate face neighbors
            //         imm_f2_normal = forward_solver->hullGeometry->faceNormal(e.halfedge().twin().face());
            double f1_area_sign = dot(f1_normal, cross(v_normal, tmp_normal)) >= 0 ? 1. : -1.; // f1 on rhs of bndN->vN
            if (forward_solver->vertex_is_stabilizable[v]){
                printf("-SV-\n");
                flow_back_boundary_on_edge_for_autodiff(bnd_normal, Edge(), v, f1_area_sign, var_positions, var_G);
            }
            else {
                for (Edge neigh_e: v.adjacentEdges()){
                    if (neigh_e != e){
                        // face_attraction_boundary[bnd_normal->f1].push_back(bnd_normal);
                        // face_attraction_boundary[bnd_normal->f2].push_back(bnd_normal);
                        printf("\n-Se %d-", neigh_e.getIndex());
                        flow_back_boundary_on_edge_for_autodiff(bnd_normal, neigh_e, v, f1_area_sign, var_positions, var_G);
                    }
                }
            }
        }
    }
    double total_area = 0.;
    autodiff::var total_area_ad = 0.;
    for (Face f: forward_solver->hullMesh->faces()){
        if (forward_solver->face_last_face[f] == f){
            total_area += face_region_area[f];
            // total_area_ad += face_region_area_ad[f];
            // std::cout << "face " << f.getIndex() << " area: " << face_region_area[f] << ",  " << face_region_area_ad[f] << std::endl;
        }
    }
    std::cout << "total face areas: " << total_area << "---" << total_area_ad << std::endl;
    
    
    // -------- gradients --------
    // if (generate_gradients){
    //     // df/dv
    //     for (Face f: forward_solver->hullMesh->faces()){
    //         // autodiff::MatrixXvar dfdv_mat;
    //         Eigen::MatrixXd dfdv_mat;
    //         if (forward_solver->face_last_face[f] == f){
    //             autodiff::VectorXvar poses_ad_vec = autodiff::VectorXvar{var_positions.reshaped()};
    //             // autodiff::VectorXvar dfdv = autodiff::gradient(face_region_area_ad[f], poses_ad_vec);
    //             Eigen::VectorXd dfdv = autodiff::gradient(face_region_area_ad[f], poses_ad_vec);
    //             dfdv_mat = dfdv.reshaped(var_positions.rows(), var_positions.cols());
    //             // df_dv_grads_ad[f] = dfdv_mat.cast<double>();
    //             df_dv_grads_ad[f] = dfdv_mat;
    //             // df/dG
    //             df_dG_grads[f] = autodiff::gradient(face_region_area_ad[f], var_G).cast<double>();
    //         }
    //     }
    // }
}


bool BoundaryBuilder::flow_back_boundary_on_edge_for_autodiff(BoundaryNormal* bnd_normal, Edge src_e, Vertex common_vertex,
                                                        double f1_area_sign, autodiff::MatrixX3var &var_positions, autodiff::Vector3var &var_G){
    // bnd_normal has to be in boundary normals of dest_e; won't assert tho for better performance
    Vertex v = common_vertex; // given as argument for better performance
    std::cout<<" . "<<std::flush;
    
    BoundaryNormal* next_bnd_normal = nullptr;
    if (forward_solver->vertex_is_stabilizable[v]){ // we are at a source
        // BoundaryNormal *curr_vertex_boundary_normal = vertex_boundary_normal[v];
        next_bnd_normal = vertex_boundary_normal[v];
        if (vertex_boundary_normal[v] == nullptr){ // first time arriving at this stable vertex; create the boundary normal
            Vector3 stable_vertex_normal = forward_solver->vertex_stable_normal[v];
            next_bnd_normal = new BoundaryNormal(stable_vertex_normal);
            // ad stuff
            next_bnd_normal->normal_ad = var_positions.row(v.getIndex()) - var_G.transpose(); // not normalizing; for more stability of gradients?
            
            vertex_boundary_normal[v] = next_bnd_normal;
            next_bnd_normal->host_v = v;
        }
        // next_bnd_normal->add_neighbor(bnd_normal);
        // bnd_normal->add_neighbor(next_bnd_normal);
        std::cout << "V!\n" << std::flush;
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
            if(e_bnd_normal.norm() == 0.){
                printf(" -src but miss- ");
                return false; // not a source for the given bnd_normal
            }
            // edge_boundary_normals[src_e].push_back(next_bnd_normal);
            edge_boundary_normals[src_e].push_back(e_bnd_normal);
            // found the source normal
            // BoundaryNormal *new_boundary_normal = new BoundaryNormal(e_bnd_normal);
            next_bnd_normal = new BoundaryNormal(e_bnd_normal);
            // ad stuff
            autodiff::Vector3var tmp_f1_normal_ad = face_normal_ad(var_positions, f1),
                                 tmp_f2_normal_ad = face_normal_ad(var_positions, f2);
            next_bnd_normal->normal_ad = intersect_arc_ray_with_arc_ad(var_positions, var_G, v, 
                                                                            bnd_normal->normal_ad,
                                                                            tmp_f1_normal_ad, tmp_f2_normal_ad, sign_change);

            // bnd_normal->add_neighbor(next_bnd_normal);
            // next_bnd_normal->add_neighbor(bnd_normal);
            next_bnd_normal->host_e = src_e;

            // region labels
            next_bnd_normal->f1 = bnd_normal->f1;
            next_bnd_normal->f2 = bnd_normal->f2;
            
            // go with the back-flow
            // printf("recursion back flow!\n");
            Vertex next_v = src_e.otherVertex(v);
            if (forward_solver->vertex_is_stabilizable[next_v]){
                std::cout << "v" << std::flush;
                flow_back_boundary_on_edge_for_autodiff(next_bnd_normal, Edge(), next_v, f1_area_sign, var_positions, var_G);
            }
            else{
                for (Edge next_src_e: next_v.adjacentEdges()){
                    if (next_src_e != src_e){
                        // non_singularity and divisive-ness will be checked inside the function
                        std::cout << "e " << next_src_e.getIndex() << std::flush;
                        bool res = flow_back_boundary_on_edge_for_autodiff(next_bnd_normal, next_src_e, next_v, f1_area_sign, var_positions, var_G);
                        if (res) // got to maximum and returning
                            break;
                    }
                }
            }
        }
        else 
            std::cout << " not src! " << std::flush;
    }
    // printf("area computes!\n");
    // printf(" computing face areas: %d, %d\n", bnd_normal->f1.getIndex(), bnd_normal->f2.getIndex());
                        
    std::cout << " F " << std::flush;
        
    if (bnd_normal->normal.norm() == 0. || next_bnd_normal == nullptr) { // if next is not found.
        std::cout << "ret \n" << std::flush;
        return false;
    }
    // printf("fetching normals and signs %d, %d\n", bnd_normal == nullptr, next_bnd_normal==nullptr);
    bool verbose = abs(face_region_area[bnd_normal->f1] - face_region_area_ad[bnd_normal->f1]) > 1e-4;
    if (verbose){
        printf("  Verboses!\n -00 NO AD %f, %f\n", face_region_area[bnd_normal->f1], face_region_area[bnd_normal->f2]);
        std::cout << "  -00    AD " << face_region_area_ad[bnd_normal->f1] << ", " << face_region_area_ad[bnd_normal->f2] << "\n";
    }
    Vector3 tmp_normal = bnd_normal->normal;
    Vector3 next_normal = next_bnd_normal->normal;
    Vector3 f1_normal = forward_solver->hullGeometry->faceNormal(bnd_normal->f1), // f1,f2 are the same along the current path to maximum
            f2_normal = forward_solver->hullGeometry->faceNormal(bnd_normal->f2);
    double curr_f1_alignment = dot(f1_normal, cross(next_normal, tmp_normal)) >= 0 ? 1. : -1.; // checking alignment again since it could change along the way
    double curr_f2_alignment = dot(f2_normal, cross(next_normal, tmp_normal)) >= 0 ? 1. : -1.;
    double f1_sign_change = f1_area_sign == curr_f1_alignment ? 1. : -1;
    double f2_sign_change = (-f1_area_sign == curr_f2_alignment) ? 1.: -1;
    autodiff::Vector3var f1_normal_ad = face_normal_ad(var_positions, bnd_normal->f1),
                         f2_normal_ad = face_normal_ad(var_positions, bnd_normal->f2);
    double f1_side_patch = triangle_patch_area_on_sphere(f1_normal, tmp_normal, next_normal),
           f2_side_patch = triangle_patch_area_on_sphere(f2_normal, tmp_normal, next_normal);
    autodiff::var f1_side_patch_ad = triangle_patch_area_on_sphere_ad(f1_normal_ad, bnd_normal->normal_ad, next_bnd_normal->normal_ad),
                  f2_side_patch_ad = triangle_patch_area_on_sphere_ad(f2_normal_ad, bnd_normal->normal_ad, next_bnd_normal->normal_ad);
    // printf("adding up patch areas\n");
    face_region_area[bnd_normal->f1] += f1_sign_change * f1_side_patch;
    face_region_area[bnd_normal->f2] += f2_sign_change * f2_side_patch;
    
    autodiff::VectorXvar poses_ad_vec = autodiff::VectorXvar{var_positions.reshaped()};
    // f1
    if (f1_side_patch != 0.){
        Eigen::VectorXd df1dv = autodiff::gradient(f1_side_patch_ad, poses_ad_vec);
        Eigen::MatrixXd df1dv_mat = df1dv.reshaped(var_positions.rows(), var_positions.cols());
        df_dv_grads_ad[bnd_normal->f1] += df1dv_mat;
        // df/dG
        df_dG_grads[bnd_normal->f1] += autodiff::gradient(f1_side_patch_ad, var_G);
    }
    // f2
    if (f2_side_patch != 0.){
        Eigen::VectorXd df2dv = autodiff::gradient(f2_side_patch_ad, poses_ad_vec);
        Eigen::MatrixXd df2dv_mat = df2dv.reshaped(var_positions.rows(), var_positions.cols());
        df_dv_grads_ad[bnd_normal->f2] += df2dv_mat;
        // df/dG
        df_dG_grads[bnd_normal->f2] += autodiff::gradient(f2_side_patch_ad, var_G);    
    }
    
    return true;

    // why is this so slow??

    if (f1_side_patch != 0.) // avoiding conditions in ad version since idk how they work yet
        face_region_area_ad[bnd_normal->f1] += f1_sign_change * f1_side_patch_ad;
    if (f2_side_patch != 0.)
        face_region_area_ad[bnd_normal->f2] += f2_sign_change * f2_side_patch_ad;
    if (verbose){
        printf("  - NO AD %f, %f\n", face_region_area[bnd_normal->f1], face_region_area[bnd_normal->f2]);
        std::cout << "  -    AD " << face_region_area_ad[bnd_normal->f1] << ", " << face_region_area_ad[bnd_normal->f2] << "\n";
        std::cout << "  - patches   " << f1_side_patch << ", " << f2_side_patch << "\n";
        std::cout << "  - AD patchs " << f1_side_patch_ad << ", " << f2_side_patch_ad << "\n";
        printf("--------------------------------------------\n");
    }
    // printf("area compute done!\n");
    return true;
    // TODO: take care of when f1,f2 on the same side when starting from saddle 
}


// gradients
void compute_df_dv_grads_autodiff(){
    FaceData<autodiff::MatrixX3var> face_region_areas_gradient_ad;

}
