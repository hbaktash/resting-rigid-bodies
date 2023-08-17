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
    // current SudoEdge will be the first, by contract
    assert(new_normal.norm() == 1.);
    if (this == next_sudo_face){
        printf("This is a terminal SudoFace");
        return nullptr;
    }
    // SudoEdge is well-defined 
    // TODO: handle new normal not being new!
    SudoFace *new_sudo_face = new SudoFace(this->host_he, new_normal, this->next_sudo_face, this);
    this->next_sudo_face->prev_sudo_face = new_sudo_face;
    this->next_sudo_face = new_sudo_face;
    // leave the twin side to the caller; to handle singular edges and to pass the twin pointers.
    return new_sudo_face;
}


// flow sf to sf and split the dest if needed, recursive and should return dest_sf2 in the end
// TODO: make it return void?? and just assert output to be dest_sf2
void RollingMarkovModel::flow_sf_to_sf(SudoFace* src_sf1, SudoFace* src_sf2, SudoFace* dest_sf1, SudoFace* dest_sf2){
    Vertex host_v = src_sf1->host_he.tipVertex(); // since src_sf is the source; already determined in initialization;
    Vector3 source_normal = vertex_stable_normal[host_v];

    Vector3 tail_intersection_normal_src = intersect_arc_ray_with_arc(source_normal, src_sf1->normal, dest_sf1->normal, dest_sf2->normal);
    if (tail_intersection_normal_src.norm() == 0.){ // no intersection
        Vector3 tail_intersection_normal_dest = intersect_arc_ray_with_arc(source_normal, dest_sf1->normal, src_sf1->normal, src_sf2->normal);
            if (tail_intersection_normal_dest.norm() == 0.){ // then no intersection on 
                return;
            }
            else { // TODO: handle and calculate the probabilities here: sf to sf; with and without splits
                Vector3 tip_intersection_normal_src = intersect_arc_ray_with_arc(source_normal, src_sf2->normal, dest_sf1->normal, dest_sf2->normal);
                if (tip_intersection_normal_src.norm() == 0.){ // then src covers all of dest; TODO calc prob of 
                    // TODO: find the probability here
                    return; // Could possibly speed up the next search here by returning the last normal
                }
                else { // need to split dest sf
                    SudoFace* new_dest_sf = dest_sf1->split_sudo_edge(tip_intersection_normal_src);
                    // TODO: find the probability here
                    flow_sf_to_sf(src_sf1, src_sf2, new_dest_sf, dest_sf2);
                }
            }
    }
    else {

    }
}

// handle single source SudoEdge
void RollingMarkovModel::outflow_sudoEdge(SudoFace* tail_sf){
    Halfedge host_he = tail_sf->host_he;
    Vertex curr_v = host_he.tipVertex(); // flow is along the halfedge
    
    // greedy approach first
    for (Halfedge he: curr_v.outgoingHalfedges()){
        SudoFace* tmp_sf = root_sudo_face[he];
        if (tmp_sf != nullptr){ // this HalfEdge has outflow from the current vertex 
            while (tmp_sf->next_sudo_face != tmp_sf){

                tmp_sf
            }
        }
    }
}

// handle single source HalfEdge
void RollingMarkovModel::outflow_halfedge(Halfedge he){
    if (edge_is_singular[he.edge()]) // if singular, won't need further outflow tracing (is already chopped by prev source halfEdges)
        return;                      // TODO: do something else or compute probabilities later?
    SudoFace* curr_sf = root_sudo_face[he];
    // SudoFace* search_init_sf = root_sudo_face[he]; // TODO: need to do smth about this for speed-up
    while (curr_sf->next_sudo_face != curr_sf){
        outflow_sudoEdge(curr_sf);
        curr_sf = curr_sf->next_sudo_face;
    }
}

// BFS on vertices/halfedges
void RollingMarkovModel::split_chain_edges(){
    // initiate the bfs queue
    bfs_list = std::list<Halfedge>();
    for (Vertex v: mesh->vertices()){
        if (vertex_stabilizablity[v]){
            for (Halfedge he: v.outgoingHalfedges()){
                bfs_list.push_back(he);
                // vertex_has_been_in_list[v] = true;
            }
        }
    }
    // do the BFS and flow the vector field
    assert(bfs_list.size() != 0); // should have at least one source/stable vertex
    while (true){
        Halfedge curr_he = bfs_list.front();
        bfs_list.pop_front();
        // TODO: do smth about curr_he
        outflow_halfedge(curr_he);
    }
}


// whether the stable point can be reached or not
void RollingMarkovModel::compute_vertex_stabilizablity(){
    vertex_stabilizablity = VertexData<bool>(*mesh, false);
    vertex_stable_normal = VertexData<Vector3>(*mesh);
    for (Vertex v: mesh->vertices()){
        if (forward_solver->vertex_is_stablizable(v))
            vertex_stabilizablity[v] = true;
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
        
        if (next_vertex.getIndex() == INVALID_IND){ // singular edge
            edge_is_singular[e] = true;
        }
        else if (next_vertex == e.secondVertex())
            edge_roll_dir[e] = 1;
        else if (next_vertex == e.firstVertex())
            edge_roll_dir[e] = -1;
        else 
            printf(" %%% This should not happen %%% \n");
        // Vertex prev_vertex = e.otherVertex(next_vertex);
        Halfedge vec_field_aligned_he = (e.secondVertex() == next_vertex) ? e.halfedge() : e.halfedge().twin();
        // // questionable; should I make the twin side null? or populate it?
        // // Going for null approach and initiating during BFS
        // initiate_root_sudo_face(e);
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


