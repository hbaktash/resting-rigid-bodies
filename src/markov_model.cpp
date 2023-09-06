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
#include "markov_model.h"

// double EPS = 1e-8;

size_t SudoFace::counter = 0;

// trivial constructors
SudoFace::SudoFace(Halfedge host_he_, Vector3 normal_, SudoFace *next_sudo_face_, SudoFace *prev_sudo_face_)
    : index(counter++) {
    host_he = host_he_;
    normal = normal_;
    next_sudo_face = next_sudo_face_;
    prev_sudo_face = prev_sudo_face_;
}

// SudoEdge::SudoEdge(Edge host_edge_, SudoFace *first_sudo_face_, SudoFace *second_sudo_face_){
//     host_edge = host_edge_;
//     first_sudo_face = first_sudo_face_;
//     second_sudo_face = second_sudo_face_;
// }

RollingMarkovModel::RollingMarkovModel(Forward3DSolver *forward_solver_){
    forward_solver = forward_solver_;
    G = forward_solver->G;
    mesh = forward_solver->hullMesh;
    geometry = forward_solver->hullGeometry;
}

RollingMarkovModel::RollingMarkovModel(ManifoldSurfaceMesh* mesh_, VertexPositionGeometry* geometry_, Vector3 G_){
    G = G_;
    mesh = mesh_;
    geometry = geometry_;
    forward_solver = new Forward3DSolver(mesh, geometry, G);
}


// split the SudoEdge starting with this SudoFace 
// TODO: source/sink assignment!
SudoFace* SudoFace::split_sudo_edge(Vector3 new_normal){
    assert(abs(new_normal.norm() - 1.) <= EPS);
    // cases that no split is required
    if (this == next_sudo_face){
        printf("This is a terminal SudoFace");
        return this;
    }
    if ((new_normal - this->normal).norm() <= EPS)
        return this;
    if ((new_normal - next_sudo_face->normal).norm() <= EPS)
        return next_sudo_face;
    // proper split is required 
    SudoFace *new_sudo_face = new SudoFace(this->host_he, new_normal, this->next_sudo_face, this);
    this->next_sudo_face->prev_sudo_face = new_sudo_face;
    this->next_sudo_face = new_sudo_face;
    return new_sudo_face;
}


// flow sf to sf and split the dest if needed
void RollingMarkovModel::flow_sf_to_sf(SudoFace* src_sf1, SudoFace* dest_sf1){
    SudoFace *src_sf2 = src_sf1->next_sudo_face,
             *dest_sf2= dest_sf1->next_sudo_face; 
    Vertex host_v = src_sf1->host_he.tipVertex(); // since src_sf is the source; already determined in initialization;
    Vector3 source_normal = vertex_stable_normal[host_v];

    // find all hit locations of the vector field
    Vector3 tail_intersection_normal_src = intersect_arc_ray_with_arc(source_normal, src_sf1->normal, dest_sf1->normal, dest_sf2->normal);
    Vector3 tip_intersection_normal_src = intersect_arc_ray_with_arc(source_normal, src_sf2->normal, dest_sf1->normal, dest_sf2->normal);
    Vector3 tail_intersection_normal_dest = intersect_arc_ray_with_arc(source_normal, dest_sf1->normal, src_sf1->normal, src_sf2->normal);
    Vector3 tip_intersection_normal_dest = intersect_arc_ray_with_arc(source_normal, dest_sf2->normal, src_sf1->normal, src_sf2->normal);
    // std::cout << "src1: (" << src_sf1->normal << ")  src2: ("<< src_sf2->normal << ") dest1: ("<< dest_sf1->normal << ") dest2: (" << dest_sf2->normal <<")\n";
    // check whether hits are inside arc segments or not
    bool tail_src_hits  = tail_intersection_normal_src.norm() >= EPS,
         tip_src_hits   = tip_intersection_normal_src.norm()  >= EPS,
         tail_dest_hits = tail_intersection_normal_dest.norm()>= EPS,
         tip_dest_hits  = tip_intersection_normal_dest.norm() >= EPS;
    
    Vertex v = dest_sf1->host_he.tailVertex();
    double total_src_angle = angle(src_sf1->normal, src_sf2->normal);
    double vertex_patch_area = vertex_gaussian_curvature[v]; // already computed
    printf("");
    // TODO: handle tail to tail cases; is it already handled???
    //       check splits and check probability pairs being generated
    // verdict: - no redundant sudoFace is being generated
    //          - redundant prob_pairs might be generated, with close to zero probability assigned to them
    if (tail_src_hits && !tip_src_hits){ // one hit. src tail
        printf("        - one hit. src tail\n");
        SudoFace* new_dest_sf = dest_sf1->split_sudo_edge(tail_intersection_normal_src);
        new_dest_sf->source_sudo_face = src_sf1;
        // assigning probability
        if (tip_dest_hits && 
            (!tail_dest_hits || (tail_intersection_normal_dest - src_sf1->normal).norm() <= EPS)){ // either no hit on the other end, or is the same as the src tail
            
            double dest_portion_angle = angle(tip_intersection_normal_dest, src_sf1->normal);
            double prob = dest_portion_angle/total_src_angle;
            sf_sf_pairs.push_back({src_sf1, new_dest_sf});
            sf_sf_probs.push_back(prob);

            // vertex to sf
            double vertex_sf_prob = patch_area(tip_intersection_normal_dest, src_sf1->normal, new_dest_sf->normal, new_dest_sf->next_sudo_face->normal);
            vertex_sf_pairs.push_back({v, new_dest_sf});
            vertex_sf_probs.push_back(vertex_sf_prob/vertex_patch_area);

        }
        else if (tail_dest_hits && 
            (!tip_dest_hits || (tip_intersection_normal_dest - src_sf1->normal).norm() <= EPS)){ // miss-aligned orientation of src and dest
            double dest_portion_angle = angle(tail_intersection_normal_dest, src_sf1->normal);
            double prob = dest_portion_angle/total_src_angle;
            sf_sf_pairs.push_back({src_sf1, dest_sf1});
            sf_sf_probs.push_back(prob);

            // vertex to sf
            double vertex_sf_prob = patch_area(tail_intersection_normal_dest, src_sf1->normal, dest_sf1->normal, tail_intersection_normal_src);
            vertex_sf_pairs.push_back({v, dest_sf1});
            vertex_sf_probs.push_back(vertex_sf_prob/vertex_patch_area);
        }
        else 
            throw std::logic_error("tip/tail intersections have gone wrong\n");
        return;
    }
    else if (!tail_src_hits && tip_src_hits){ // one hit. src tip
        printf("        - one hit. src tip\n");
        SudoFace* new_dest_sf = dest_sf1->split_sudo_edge(tip_intersection_normal_src);
        new_dest_sf->source_sudo_face = src_sf2;
        // assigning probability
        // assert((tail_dest_hits && !tip_dest_hits) || (!tail_dest_hits && tip_dest_hits)); // exactly one should hit
        if (tip_dest_hits && 
            (!tail_dest_hits || (tail_intersection_normal_dest - src_sf2->normal).norm() <= EPS)){ // aligned orientation of src and dest
            double dest_portion_angle = angle(tip_intersection_normal_dest, src_sf2->normal);
            double prob = dest_portion_angle/total_src_angle;
            sf_sf_pairs.push_back({src_sf1, new_dest_sf});
            sf_sf_probs.push_back(prob);

            // vertex to sf
            double vertex_sf_prob = patch_area(tip_intersection_normal_dest, src_sf2->normal, new_dest_sf->normal, new_dest_sf->next_sudo_face->normal);
            vertex_sf_pairs.push_back({v, new_dest_sf});
            vertex_sf_probs.push_back(vertex_sf_prob/vertex_patch_area);
        }
        else if (tail_dest_hits && 
            (!tip_dest_hits || (tip_intersection_normal_dest - src_sf2->normal).norm() <= EPS)){ // miss-aligned orientation of src and dest
            double dest_portion_angle = angle(tail_intersection_normal_dest, src_sf2->normal);
            double prob = dest_portion_angle/total_src_angle;
            sf_sf_pairs.push_back({src_sf1, dest_sf1});
            sf_sf_probs.push_back(prob);
            
            // vertex to sf
            double tmp_patch_area = patch_area(tail_intersection_normal_dest, src_sf2->normal, dest_sf1->normal, tip_intersection_normal_src);
            vertex_sf_pairs.push_back({v, dest_sf1});
            vertex_sf_probs.push_back(tmp_patch_area/vertex_patch_area);

            printf("            (v,sf): %d, %d   prob: %f / %f\n", v.getIndex(), dest_sf1->index, tmp_patch_area, vertex_patch_area);
        }
        else 
            throw std::logic_error("tip/tail intersections have gone wrong\n");
        return;
    }
    else if (!tail_src_hits && !tip_src_hits){ //  no hits!
        printf("        - not hits\n");
        if (tail_dest_hits && tip_dest_hits){ // all of dest sudoEdge is covered; one could hit because of pre/next SudoEdge overlaps
            // assigning probability
            double dest_portion_angle = angle(tip_intersection_normal_dest, tail_intersection_normal_dest);
            double prob = dest_portion_angle/total_src_angle;
            sf_sf_pairs.push_back({src_sf1, dest_sf1});
            sf_sf_probs.push_back(prob);
            
            // vertex to sf
            double vertex_sf_prob = patch_area(tip_intersection_normal_dest, tail_intersection_normal_dest, dest_sf1->normal, dest_sf1->next_sudo_face->normal);
            vertex_sf_pairs.push_back({v, dest_sf1});
            vertex_sf_probs.push_back(vertex_sf_prob/vertex_patch_area);
            return;
        }
        else { // total miss on dest sudoEdge; 
            // assigning probability
            // is zero by default
            return;
        }
    }
    else if (tail_src_hits && tip_src_hits){ // two hits!!
        printf("        - two hits!!\n");
        SudoFace *new_dest_sf1, *new_dest_sf2;
        // need to check alignment here
        if ((tail_intersection_normal_src - dest_sf1->normal).norm() <= (tip_intersection_normal_src - dest_sf1->normal).norm()){ // src and dest have aligned orientations 
            new_dest_sf1 = dest_sf1->split_sudo_edge(tail_intersection_normal_src);
            new_dest_sf2 = new_dest_sf1->split_sudo_edge(tip_intersection_normal_src);
            new_dest_sf1->source_sudo_face = src_sf1;
            new_dest_sf2->source_sudo_face = src_sf2;
        }
        else {
            new_dest_sf1 = dest_sf1->split_sudo_edge(tip_intersection_normal_src);
            new_dest_sf2 = new_dest_sf1->split_sudo_edge(tail_intersection_normal_src);
            new_dest_sf1->source_sudo_face = src_sf2;
            new_dest_sf2->source_sudo_face = src_sf1;
        }
        // assigning probabilities
        double prob = 1.;
        sf_sf_pairs.push_back({src_sf1, new_dest_sf1});
        sf_sf_probs.push_back(prob);
        
        // vertex to sf
        double tmp_patch_area = patch_area(src_sf1->normal, src_sf2->normal, 
                                           new_dest_sf1->normal, new_dest_sf1->next_sudo_face->normal);
        vertex_sf_pairs.push_back({v, new_dest_sf1});
        vertex_sf_probs.push_back(tmp_patch_area/vertex_patch_area);
        printf("            (v,sf): %d, %d   prob: %f / %f\n", v.getIndex(), new_dest_sf1->index, tmp_patch_area,vertex_patch_area);
        return;
    }
    else printf("        - WTF hits??? \n");; //shouldnt get here!
}


void RollingMarkovModel::flow_he_to_he(Halfedge src, Halfedge dest){
    SudoFace *curr_src_sf = root_sudo_face[src];
    SudoFace *root_dest_sf = root_sudo_face[dest]; // never changes; even due to splits;
    assert(root_dest_sf != nullptr && curr_src_sf != nullptr);
    int src_cnt = 0, dest_cnt = 0;
    while (curr_src_sf->next_sudo_face != curr_src_sf) {
        SudoFace *curr_dest_sf = root_dest_sf;
        while(curr_dest_sf->next_sudo_face != curr_dest_sf){
            printf("     ** flow src sf %d, to dest %d \n", curr_src_sf->index, curr_dest_sf->index);
            flow_sf_to_sf(curr_src_sf, curr_dest_sf);
            curr_dest_sf = curr_dest_sf->next_sudo_face;
            dest_cnt++;
        }
        curr_src_sf = curr_src_sf->next_sudo_face;
        src_cnt++;
    }
}

// handle single sink HalfEdge
void RollingMarkovModel::process_halfedge(Halfedge he){
    if (he_processed[he])
        return;
    if (vertex_is_stabilizable[he.tailVertex()]) // No split needed. The edge is fully reachable from the stable tail vertex
        return;
    // else
        // assert(forward_solver->next_rolling_vertex(he.edge()) == he.tipVertex()); // the HalfEdge aligns with the flow
    printf("  -- processing he %d, %d  f, tf: %d,%d\n", he.tailVertex().getIndex(), he.tipVertex().getIndex(),
                                                        he.face().getIndex(), he.twin().face().getIndex());
    
    Vertex v = he.tailVertex();
    for (Halfedge src_he: v.incomingHalfedges()){
        // printf("    & checking possible src_he %d, %d,   singular: %d, has root sf %d \n", src_he.tailVertex().getIndex(), src_he.tipVertex().getIndex(), edge_is_singular[src_he.edge()], root_sudo_face[src_he] != nullptr);
        if (root_sudo_face[src_he] != nullptr && !edge_is_singular[src_he.edge()]){ // src_he is a source in this vertex
            process_halfedge(src_he);
            printf("    * flowing from %d,%d  f, tf: %d,%d \n", src_he.tailVertex().getIndex(), 
                                                                 src_he.tipVertex().getIndex(),
                                                                 src_he.face().getIndex(),
                                                                 src_he.twin().face().getIndex());
            flow_he_to_he(src_he, he);
        }
    }
    he_processed[he] = true;
}


void RollingMarkovModel::empty_prob_vectors(){
    vertex_sf_pairs.clear();
    vertex_sf_probs.clear();
    sf_sf_pairs.clear();
    sf_sf_probs.clear();
    sf_face_pairs.clear();
    sf_face_probs.clear();
}


// recursion starting from singular/stable edges
void RollingMarkovModel::split_chain_edges_and_build_probability_pairs(){
    // clean up probability pair vectors
    empty_prob_vectors();

    // DP to avoid spliting a HalfEdge twice
    he_processed = HalfedgeData<bool>(*mesh, false);
    // use singular edges as starting seeds the recursion
    printf("finding terminal seeds \n");
    for (Edge e: mesh->edges()){ 
        if (edge_is_singular[e]){ // go back up from singular edges; split any edge if needed
            printf("  -- at singular edge %d, %d \n", e.firstVertex().getIndex(), e.secondVertex().getIndex());
            Vertex v1 = e.firstVertex(), 
                   v2 = e.secondVertex();
            Halfedge he = e.halfedge(),
                     he_twin = e.halfedge().twin();
            SudoFace *sf = root_sudo_face[he],
                     *sf_twin = root_sudo_face[he_twin]; 
            if (!vertex_is_stabilizable[v1])  // v1 is not a source/stable
                termilar_hes.push_back(he);
            else {
                he_processed[he] = true;
                
                // assigning probabilities
                // vertex -> edge
                Vector3 v1_stable_normal = vertex_stable_normal[v1];
                double v1f1f2_patch_area = triangle_patch_area_on_sphere(v1_stable_normal, 
                                                                         geometry->faceNormal(he.face()),
                                                                         geometry->faceNormal(he.twin().face()));
                vertex_sf_pairs.push_back({v1, sf});
                vertex_sf_probs.push_back(v1f1f2_patch_area/vertex_gaussian_curvature[v1]);
                // printf("added v_sf pair: %d  %d,%d\n", v1.getIndex(), sf->host_he.tailVertex().getIndex(),
                //                                                       sf->host_he.tipVertex().getIndex());
                // printf("    patch area, vertex G_Curvature: %f , %f\n", v1f1f2_patch_area, vertex_gaussian_curvature[v1]);
            }
            
            if (!vertex_is_stabilizable[v2]) // v2 is not a source/stable
                termilar_hes.push_back(he_twin);
            else {

                he_processed[he_twin] = true;
                // vertex -> edge
                Vector3 v2_stable_normal = vertex_stable_normal[v2];
                double v2f1f2_patch_area = triangle_patch_area_on_sphere(v2_stable_normal, 
                                                                         geometry->faceNormal(he.face()),
                                                                         geometry->faceNormal(he_twin.face()));
                vertex_sf_pairs.push_back({v2, sf_twin});
                vertex_sf_probs.push_back(v2f1f2_patch_area/vertex_gaussian_curvature[v2]);
            }
        }
    }
    // Recursive DFS from every starting seed edge (singular edge)
    printf("processing halfedges\n");
    for (Halfedge he: termilar_hes){
        process_halfedge(he);
    }

    build_sf_face_pairs();
}

void RollingMarkovModel::build_sf_face_pairs(){
    for (Halfedge he: mesh->halfedges()){
        if (edge_is_singular[he.edge()]){
            Face f1 = he.face(),
                 f2 = he.twin().face();
            Vector3 stable_normal = edge_stable_normal[he.edge()];
            SudoFace *curr_sf = root_sudo_face[he];
            if (curr_sf != nullptr){
                while (curr_sf->next_sudo_face != curr_sf) {
                    Vector3 curr_normal1 = curr_sf->normal,
                            curr_normal2 = curr_sf->next_sudo_face->normal;
                    double f1_prob = arc_portion(stable_normal, curr_normal1, curr_normal2);
                    sf_face_pairs.push_back({curr_sf, f1});
                    sf_face_probs.push_back(f1_prob);
                    sf_face_pairs.push_back({curr_sf, f2});
                    sf_face_probs.push_back(1. - f1_prob);
                    curr_sf = curr_sf->next_sudo_face;
                }
            }
        }
    }
}

// whether the stable point can be reached or not
void RollingMarkovModel::compute_vertex_stabilizablity(){
    vertex_is_stabilizable = VertexData<bool>(*mesh, false);
    vertex_stable_normal = VertexData<Vector3>(*mesh);
    for (Vertex v: mesh->vertices()){
        if (forward_solver->vertex_is_stablizable(v))
            vertex_is_stabilizable[v] = true;
        vertex_stable_normal[v] = (geometry->inputVertexPositions[v] - G).normalize();
    }
}

// initialize 
void RollingMarkovModel::initiate_root_sudo_face(Halfedge he){
    // Halfedge he = e.halfedge();
    Vertex v1 = he.tailVertex(),
           v2 = he.tipVertex();
    SudoFace *sf1 = new SudoFace(he, geometry->faceNormal(he.face()), nullptr, nullptr),
             *sf2 = new SudoFace(he, geometry->faceNormal(he.twin().face()), nullptr, nullptr);
    sf1->next_sudo_face = sf2;
    sf1->prev_sudo_face = sf1;
    sf2->next_sudo_face = sf2;
    sf2->prev_sudo_face = sf1;

    root_sudo_face[he] = sf1;
}

// edge rolls to a face if singular; else rolls to a vertex 
void RollingMarkovModel::compute_edge_singularity_and_init_source_dir(){
    edge_is_singular = EdgeData<bool>(*mesh, false); // ~ is stable, by fwdSolver function names
    // initial assignment of real faces as SudoFaces
    root_sudo_face = HalfedgeData<SudoFace*>(*mesh, nullptr);
    for (Edge e: mesh->edges()){
        Vertex next_vertex = forward_solver->next_rolling_vertex(e);
        
        if (next_vertex.getIndex() == INVALID_IND){ // singular edge
            edge_is_singular[e] = true;
            initiate_root_sudo_face(e.halfedge());
            initiate_root_sudo_face(e.halfedge().twin());
        }
        else {
            Halfedge vec_field_aligned_he = (e.firstVertex() == next_vertex) ? e.halfedge().twin() : e.halfedge();
            // so singular he is also aligned with the vector field
            initiate_root_sudo_face(vec_field_aligned_he);
            // printf(" inited sf for non singular he %d, %d\n", vec_field_aligned_he.tailVertex().getIndex(), vec_field_aligned_he.tipVertex().getIndex());
        }
    }
}

// is 0 if the normal is unreachable; and if non-singular: the normal doesnt fall on the edge (edge too short)
void RollingMarkovModel::compute_edge_stable_normals(){
    // edge_is_singular should be populated and initiated before
    if (edge_is_singular.size() == 0) throw std::logic_error("edge_is_singular should be called before this.\n");
    
    Vector3 zero_vec({0.,0.,0.});
    edge_stable_normal = EdgeData<Vector3>(*mesh, zero_vec);
    for (Edge e: mesh->edges()){
        if (edge_is_singular[e]){ // not cheking stabilizability ~= reachablility
            Vector3 A = geometry->inputVertexPositions[e.firstVertex()],
                    B = geometry->inputVertexPositions[e.secondVertex()];
            edge_stable_normal[e] = point_to_segment_normal(G, A, B).normalize();
        }
    }
}

// just call all the pre-compute initializations
void RollingMarkovModel::initialize_pre_computes(){
    SudoFace::counter = 0;
    compute_vertex_stabilizablity();
    compute_vertex_gaussian_curvatures();
    compute_edge_singularity_and_init_source_dir();
    compute_edge_stable_normals();
}

// pre-compute vertex gaussian curvatures
void RollingMarkovModel::compute_vertex_gaussian_curvatures(){
    vertex_gaussian_curvature = VertexData<double>(*mesh, 0.);
    for (Vertex v: mesh->vertices()){
        vertex_gaussian_curvature[v] = gaussian_curvature(v, *geometry);
        // printf(" Gaussian Curvature at %d, is %f\n", v.getIndex(), vertex_gaussian_curvature[v]);
    }
}

//
double RollingMarkovModel::get_he_face_probability(Halfedge he){
    if (edge_is_singular[he.edge()]){
        Vector3 f1_normal = geometry->faceNormal(he.face()),
                f2_normal = geometry->faceNormal(he.twin().face()),
                stable_normal = edge_stable_normal[he.edge()];
        assert(stable_normal.norm() >= EPS); // since edge is singular
        return arc_portion(stable_normal, f1_normal, f2_normal);
    }
    else 
        return 0.;
}

void RollingMarkovModel::print_prob_pairs(){
    printf("--printing probability pairs--\n\n");
    printf("Vertex - SF:he -  probs: \n");
    for (int i = 0; i < vertex_sf_pairs.size(); i++){
        std::pair<Vertex, SudoFace*> v_sf_pair = vertex_sf_pairs[i];
        Vertex v = v_sf_pair.first;
        SudoFace* sf = v_sf_pair.second;
        double v_sf_prob = vertex_sf_probs[i];
        printf("  %d  -  (%d) %d,%d  -  %f  \n", v.getIndex(), sf->index,
                                            sf->host_he.tailVertex().getIndex(), 
                                            sf->host_he.tipVertex().getIndex(), 
                                            v_sf_prob);
    }
    printf("\nSF1:he   -   SF2:he   -   probs: \n");
    for (int i = 0; i < sf_sf_pairs.size(); i++){
        std::pair<SudoFace*, SudoFace*> sf_sf_pair = sf_sf_pairs[i];
        SudoFace *sf1 = sf_sf_pair.first,
                 *sf2 = sf_sf_pair.second;
        double sf_sf_prob = sf_sf_probs[i];
        printf("  (%d) %d,%d  ->  (%d) %d,%d  =  %f  \n", sf1->index, 
                                            sf1->host_he.tailVertex().getIndex(), 
                                            sf1->host_he.tipVertex().getIndex(), sf2->index,
                                            sf2->host_he.tailVertex().getIndex(), 
                                            sf2->host_he.tipVertex().getIndex(), 
                                            sf_sf_prob);
    }
    printf("\nSF: he   -   face   -   probs: \n");
    for (int i = 0; i < sf_face_pairs.size(); i++){
        std::pair<SudoFace*, Face> sf_face_pair = sf_face_pairs[i];
        SudoFace *sf = sf_face_pair.first;
        Face f = sf_face_pair.second;
        double sf_face_prob = sf_face_probs[i];
        printf("  sf(%d) %d,%d  ->   f(%d)  =  %f  \n", sf->index, 
                                            sf->host_he.tailVertex().getIndex(), 
                                            sf->host_he.tipVertex().getIndex(), 
                                            f.getIndex(), 
                                            sf_face_prob);
    }
}