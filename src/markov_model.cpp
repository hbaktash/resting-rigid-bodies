#include "markov_model.h"



// trivial constructors
SudoFace::SudoFace(Halfedge host_he_, Vector3 normal_, SudoFace *next_sudo_face_, SudoFace *prev_sudo_face_){
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
// TODO: decide on twin he assignment; null/potent for sink side
SudoFace* SudoFace::split_sudo_edge(Vector3 new_normal){
    assert(new_normal.norm() == 1.);
    // cases that no split is required
    if (this == next_sudo_face){
        printf("This is a terminal SudoFace");
        return this;
    }
    if (new_normal == this->normal)
        return this;
    if (new_normal == next_sudo_face->normal)
        return next_sudo_face;
    // proper split is required 
    SudoFace *new_sudo_face = new SudoFace(this->host_he, new_normal, this->next_sudo_face, this);
    this->next_sudo_face->prev_sudo_face = new_sudo_face;
    this->next_sudo_face = new_sudo_face;
    return new_sudo_face;
}


// flow sf to sf and split the dest if needed, recursive and should return dest_sf2 in the end
// TODO: make it return void?? and just assert output to be dest_sf2
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
    // check whether hits are inside arc segments or not
    bool tail_src_hits  = tail_intersection_normal_src.norm() != 0.,
         tip_src_hits   = tip_intersection_normal_src.norm() != 0.,
         tail_dest_hits = tail_intersection_normal_dest.norm() != 0.,
         tip_dest_hits  = tip_intersection_normal_dest.norm() != 0.;
    if (tail_src_hits && !tip_src_hits){ // one hit. src tail
        SudoFace* new_dest_sf = dest_sf1->split_sudo_edge(tail_intersection_normal_src);
        new_dest_sf->source_sudo_face = src_sf1;
        // TODO: check dest hits to get probability
        return;
    }
    else if (!tail_src_hits && tip_src_hits){ // one hit. src tip
        SudoFace* new_dest_sf = dest_sf1->split_sudo_edge(tip_intersection_normal_src);
        new_dest_sf->source_sudo_face = src_sf2;
        // TODO: check dest hits to get probability
        return;
    }
    else if (!tail_src_hits && !tip_src_hits){ //  no hits!
        if (tail_dest_hits){ // all of dest sudoEdge is covered 
            assert(tip_dest_hits); // either none hit or both should
            //TODO: check dest hits to get probability
            return;
        }
        else { // total miss on dest sudoEdge
            assert(!tip_dest_hits); // either none hit or both should
            //TODO: probability is zero, assign it smwhere
            return;
        }
    }
    else if (tail_src_hits && tip_src_hits){ // two hits!!
        SudoFace *new_dest_sf1, *new_dest_sf2;
        // need to check alignment here
        if ((tail_intersection_normal_src - dest_sf1->normal).norm() <= (tip_intersection_normal_src - dest_sf1->normal).norm()){ // tail-tip assignments agree on src and dest 
            new_dest_sf1 = dest_sf1->split_sudo_edge(tail_intersection_normal_src);
            new_dest_sf2 = new_dest_sf1->split_sudo_edge(tip_intersection_normal_src);
            new_dest_sf1->source_sudo_face = src_sf1;
            new_dest_sf2->source_sudo_face = src_sf2;
            //TODO: assign the probability here
        }
        else {
            new_dest_sf1 = dest_sf1->split_sudo_edge(tip_intersection_normal_src);
            new_dest_sf2 = new_dest_sf1->split_sudo_edge(tail_intersection_normal_src);
            new_dest_sf1->source_sudo_face = src_sf2;
            new_dest_sf2->source_sudo_face = src_sf1;
            //TODO: assign the probability here
        }
        return;
    }
    else; //shouldnt get here!
}

void RollingMarkovModel::flow_he_to_he(Halfedge src, Halfedge dest){
    SudoFace *curr_src_sf = root_sudo_face[src];
    SudoFace *root_dest_sf = root_sudo_face[dest]; // never changes; even due to splits;
    assert(root_dest_sf != nullptr && curr_src_sf != nullptr);
    while (curr_src_sf->next_sudo_face != curr_src_sf) {
        SudoFace *curr_dest_sf = root_dest_sf;
        while(curr_dest_sf->next_sudo_face != curr_dest_sf){
            flow_sf_to_sf(curr_src_sf, curr_dest_sf);
            curr_dest_sf = curr_dest_sf->next_sudo_face;
        }
        curr_src_sf = curr_src_sf->next_sudo_face;
    }
}

// handle single sink HalfEdge
void RollingMarkovModel::process_halfedge(Halfedge he){
    if (vertex_is_stabilizable[he.tailVertex()]) // the edge is fully reachable from the stable tail vertex
        return;
    else
        assert(forward_solver->next_rolling_vertex(he.edge()) == he.tipVertex()); // the he aligns with the flow
    if (he_processed[he])
        return;
    
    Vertex v = he.tailVertex();
    for (Halfedge src_he: v.incomingHalfedges()){
        if (root_sudo_face[src_he] != nullptr && !he_processed[src_he]){ // HalfEdge is a source in this vertex
            process_halfedge(src_he);
        }
        flow_he_to_he(src_he, he);
    }
    he_processed[he] = true;
}

// recursion starting from singular/stable edges
void RollingMarkovModel::split_chain_edges(){
    // DP to avoid spliting a HalfEdge twice
    he_processed = HalfedgeData<bool>(*mesh, false);
    // use singular edges as starting seeds the recursion
    for (Edge e: mesh->edges()){ 
        if (edge_is_singular[e]){ // go back up from singular edges; split any edge if needed
            Vertex v1 = e.firstVertex(), 
                   v2 = e.secondVertex();
            if (!vertex_is_stabilizable[v1])  // v1 is not a source/stable
                termilar_hes.push_back(e.halfedge());
            else
                he_processed[e.halfedge()] = true;
            
            if (!vertex_is_stabilizable[v2]) // v2 is not a source/stable
                termilar_hes.push_back(e.halfedge().twin());
            else
                he_processed[e.halfedge().twin()] = true;

            // TODO: find probabilities for else cases immediately
        }
    }
    // recursive DFS from every starting seed edge (singular edge)
    for (Halfedge he: termilar_hes){
        process_halfedge(he);
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
            //  *sf1_twin = new SudoFace(he.twin(), geometry->faceNormal(he.twin().face()), nullptr, nullptr),
            //  *sf2_twin = new SudoFace(he.twin(), geometry->faceNormal(he.face()), nullptr, nullptr);
    sf1->next_sudo_face = sf2;
    sf1->prev_sudo_face = sf1;
    sf2->next_sudo_face = sf2;
    sf2->prev_sudo_face = sf1;

    // sf1_twin->next_sudo_face = sf1_twin;
    // sf1_twin->prev_sudo_face = sf2_twin;
    // sf2_twin->next_sudo_face = sf1_twin;
    // sf2_twin->prev_sudo_face = sf2_twin;

    // sf1->twin = sf1_twin;
    // sf1_twin->twin = sf1;
    // sf2->twin = sf2_twin;
    // sf2_twin->twin = sf2;
}

// edge rolls to a face if singular; else rolls to a vertex 
void RollingMarkovModel::compute_edge_singularity_and_init_source_dir(){
    edge_is_singular = EdgeData<bool>(*mesh, false); // ~ is stable, by fwdSolver function names
    edge_roll_dir = EdgeData<int>(*mesh, 0);
    // initial assignment of real faces as SudoFaces
    root_sudo_face = HalfedgeData<SudoFace*>(*mesh, nullptr);
    for (Edge e: mesh->edges()){
        Vertex next_vertex = forward_solver->next_rolling_vertex(e);
        
        if (next_vertex.getIndex() == INVALID_IND) // singular edge
            edge_is_singular[e] = true;
        else if (next_vertex == e.secondVertex())
            edge_roll_dir[e] = 1;
        else if (next_vertex == e.firstVertex())
            edge_roll_dir[e] = -1;
        else 
            printf(" %%% This should not happen %%% \n");
        
        Halfedge vec_field_aligned_he = (e.firstVertex() == next_vertex) ? e.halfedge().twin() : e.halfedge();
        // so singular he is also aligned with the vector field
        initiate_root_sudo_face(vec_field_aligned_he);
    }
}

// is 0 if the normal is unreachable; and if non-singular: the normal doesnt fall on the edge (edge too short)
void RollingMarkovModel::compute_edge_stable_normal(){
    // edge_is_singular should be populated and initiated before
    if (edge_is_singular.size() == 0) throw std::logic_error("edge_is_singular should be called before this.\n");
    
    Vector3 zero_vec({0.,0.,0.});
    EdgeData<Vector3> edge_stable_normal(*mesh, zero_vec);
    for (Edge e: mesh->edges()){
        if (edge_is_singular[e] && forward_solver->edge_is_stablizable(e)){ // stabilizable ~= reachable
            Vector3 A = geometry->inputVertexPositions[e.firstVertex()],
                    B = geometry->inputVertexPositions[e.secondVertex()];
            edge_stable_normal[e] = point_to_segment_normal(G, A, B).normalize();
        }
    }
}


