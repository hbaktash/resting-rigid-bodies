#include "markov_model.h"

// double EPS = 1e-8;

size_t SudoFace::counter = 0;

// trivial constructors
SudoFace::SudoFace(Halfedge host_he_, Vector3 normal_, SudoFace *next_sudo_face_, SudoFace *prev_sudo_face_): index(counter++){
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
    double vertex_patch_area = geometry->vertexGaussianCurvature(v);

    if (tail_src_hits && !tip_src_hits){ // one hit. src tail
        SudoFace* new_dest_sf = dest_sf1->split_sudo_edge(tail_intersection_normal_src);
        new_dest_sf->source_sudo_face = src_sf1;
        // assigning probability
        assert((tail_dest_hits && !tip_dest_hits) || (!tail_dest_hits && tip_dest_hits)); // exactly one should hit
        if (tip_dest_hits){ // aligned orientation of src and dest
            double dest_portion_angle = angle(tip_intersection_normal_dest, src_sf1->normal);
            double prob = dest_portion_angle/total_src_angle;
            sf_sf_pairs.push_back({src_sf1, new_dest_sf});
            sf_sf_probs.push_back(prob);

            // vertex to sf
            double vertex_sf_prob = patch_area(tip_intersection_normal_dest, src_sf1->normal, new_dest_sf->normal, new_dest_sf->next_sudo_face->normal);
            vertex_sf_pairs.push_back({v, new_dest_sf});
            vertex_sf_probs.push_back(vertex_sf_prob/vertex_patch_area);

        }
        if (tail_dest_hits){ // miss-aligned orientation of src and dest
            double dest_portion_angle = angle(tail_intersection_normal_dest, src_sf1->normal);
            double prob = dest_portion_angle/total_src_angle;
            sf_sf_pairs.push_back({src_sf1, dest_sf1});
            sf_sf_probs.push_back(prob);

            // vertex to sf
            double vertex_sf_prob = patch_area(tail_intersection_normal_dest, src_sf1->normal, dest_sf1->normal, dest_sf1->next_sudo_face->normal);
            vertex_sf_pairs.push_back({v, dest_sf1});
            vertex_sf_probs.push_back(vertex_sf_prob/vertex_patch_area);
        }
        return;
    }
    else if (!tail_src_hits && tip_src_hits){ // one hit. src tip
        SudoFace* new_dest_sf = dest_sf1->split_sudo_edge(tip_intersection_normal_src);
        new_dest_sf->source_sudo_face = src_sf2;
        // assigning probability
        assert((tail_dest_hits && !tip_dest_hits) || (!tail_dest_hits && tip_dest_hits)); // exactly one should hit
        if (tip_dest_hits){ // aligned orientation of src and dest
            double dest_portion_angle = angle(tip_intersection_normal_dest, src_sf2->normal);
            double prob = dest_portion_angle/total_src_angle;
            sf_sf_pairs.push_back({src_sf1, new_dest_sf});
            sf_sf_probs.push_back(prob);
            // TODO vertex sf prob
        }
        if (tail_dest_hits){ // miss-aligned orientation of src and dest
            double dest_portion_angle = angle(tail_intersection_normal_dest, src_sf2->normal);
            double prob = dest_portion_angle/total_src_angle;
            sf_sf_pairs.push_back({src_sf1, dest_sf1});
            sf_sf_probs.push_back(prob);
            // TODO vertex sf prob
        }
        return;
    }
    else if (!tail_src_hits && !tip_src_hits){ //  no hits!
        if (tail_dest_hits){ // all of dest sudoEdge is covered 
            assert(tip_dest_hits); // either none hit or both should
            // assigning probability
            double dest_portion_angle = angle(tip_intersection_normal_dest, tail_intersection_normal_dest);
            double prob = dest_portion_angle/total_src_angle;
            sf_sf_pairs.push_back({src_sf1, dest_sf1});
            sf_sf_probs.push_back(prob);
            // TODO vertex sf prob
            return;
        }
        else { // total miss on dest sudoEdge
            assert(!tip_dest_hits); // either none hit or both should
            // assigning probability
            // is zero by default
            // TODO vertex sf prob
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
    if (he_processed[he])
        return;
    if (vertex_is_stabilizable[he.tailVertex()]) // No split needed. The edge is fully reachable from the stable tail vertex
        return;
    // else
        // assert(forward_solver->next_rolling_vertex(he.edge()) == he.tipVertex()); // the HalfEdge aligns with the flow
    printf("  -- processing he %d, %d\n", he.tailVertex().getIndex(), he.tipVertex().getIndex());
    
    Vertex v = he.tailVertex();
    for (Halfedge src_he: v.incomingHalfedges()){
        printf("    & checking possible src_he %d, %d,   singular: %d, has root sf %d \n", src_he.tailVertex().getIndex(), src_he.tipVertex().getIndex(), edge_is_singular[src_he.edge()], root_sudo_face[src_he] != nullptr);
        if (root_sudo_face[src_he] != nullptr && !edge_is_singular[src_he.edge()]){ // src_he is a source in this vertex
            process_halfedge(src_he);
            printf("    * flowing from %d,%d \n", src_he.tailVertex().getIndex(), src_he.tipVertex().getIndex());
            flow_he_to_he(src_he, he);
        }
    }
    he_processed[he] = true;
}

// recursion starting from singular/stable edges
void RollingMarkovModel::split_chain_edges(){
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
            if (!vertex_is_stabilizable[v1])  // v1 is not a source/stable
                termilar_hes.push_back(he);
            else {
                he_processed[he] = true;
                double face_roll_prob = get_he_face_probability(he);
                sf_face_pairs.push_back({root_sudo_face[he], he.face()});
                sf_face_probs.push_back(face_roll_prob);
                sf_face_pairs.push_back({root_sudo_face[he], he_twin.face()});
                sf_face_probs.push_back(1. - face_roll_prob);
            }
            
            if (!vertex_is_stabilizable[v2]) // v2 is not a source/stable
                termilar_hes.push_back(he_twin);
            else {
                he_processed[he_twin] = true;
                double face_roll_prob = get_he_face_probability(he_twin);
                sf_face_pairs.push_back({root_sudo_face[he_twin], he_twin.face()});
                sf_face_probs.push_back(face_roll_prob);
                sf_face_pairs.push_back({root_sudo_face[he_twin], he.face()});
                sf_face_probs.push_back(1. - face_roll_prob);
            }
            // TODO: assign vertex edge probabilities for singular edges
        }
    }
    // Recursive DFS from every starting seed edge (singular edge)
    printf("processing halfedges\n");
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
            printf(" inited sf for non singular he %d, %d\n", vec_field_aligned_he.tailVertex().getIndex(), vec_field_aligned_he.tipVertex().getIndex());
        }
    }
}

// is 0 if the normal is unreachable; and if non-singular: the normal doesnt fall on the edge (edge too short)
void RollingMarkovModel::compute_edge_stable_normals(){
    // edge_is_singular should be populated and initiated before
    if (edge_is_singular.size() == 0) throw std::logic_error("edge_is_singular should be called before this.\n");
    
    Vector3 zero_vec({0.,0.,0.});
    EdgeData<Vector3> edge_stable_normal(*mesh, zero_vec);
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
    compute_vertex_stabilizablity();
    compute_edge_singularity_and_init_source_dir();
    compute_edge_stable_normals();
}

//
double RollingMarkovModel::get_he_face_probability(Halfedge he){
    if (edge_is_singular[he.edge()]){
        Vector3 f1_normal = geometry->faceNormal(he.face()),
                f2_normal = geometry->faceNormal(he.twin().face()),
                stable_normal = edge_stable_normal[he.edge()];
        assert(stable_normal.norm() >= EPS);
        double total_angle = angle(f1_normal, f2_normal),
               sn_f1_angle  = angle(f1_normal, stable_normal),
               sn_f2_angle  = angle(f2_normal, stable_normal);
        if (sn_f1_angle <= total_angle && sn_f2_angle <= total_angle)
            return sn_f1_angle/total_angle;
        else if (sn_f1_angle >= total_angle)
            return 1.; // goes to f1
        else if (sn_f2_angle >= total_angle)
            return 0.; // goes to f2
        else
            assert(false); // not possible
    }
    else 
        return 0.;
}