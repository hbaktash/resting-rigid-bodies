#include "markov_model.h"



// trivial constructors
SudoFace::SudoFace(Edge host_edge_, Vector3 normal_, SudoFace *next_sudo_face_, SudoFace *prev_sudo_face_){
    host_edge = host_edge_;
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
    // current SudoEdge will be the first, by contract
    if (this == next_sudo_face){
        printf("This is a terminal SudoFace");
        return nullptr;
    }

    // SudoEdge is well-defined 
    SudoFace *new_sudo_face = new SudoFace(host_he, new_normal, next_sudo_face, this);
    next_sudo_face->prev_sudo_face = new_sudo_face;
    this->next_sudo_face = new_sudo_face;
    // leave the twin side to the caller; to handle singular edges and to pass the twin pointers.
    return new_sudo_face;
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
void RollingMarkovModel::initiate_root_sudo_face(Edge e){
    Halfedge he = e.halfedge();
    Vertex v1 = he.tailVertex(),
           v2 = he.tipVertex();
    SudoFace *sf1 = new SudoFace(he, geometry->faceNormal(he.face()), nullptr, nullptr),
             *sf2 = new SudoFace(he, geometry->faceNormal(he.twin().face()), nullptr, nullptr),
             *sf1_twin = new SudoFace(he.twin(), geometry->faceNormal(he.twin().face()), nullptr, nullptr),
             *sf2_twin = new SudoFace(he.twin(), geometry->faceNormal(he.face()), nullptr, nullptr);
    sf1->next_sudo_face = sf2;
    sf1->prev_sudo_face = sf1;
    sf2->next_sudo_face = sf2;
    sf2->prev_sudo_face = sf1;

    sf1_twin->next_sudo_face = sf1_twin;
    sf1_twin->prev_sudo_face = sf2_twin;
    sf2_twin->next_sudo_face = sf1_twin;
    sf2_twin->prev_sudo_face = sf2_twin;

    sf1->twin = sf1_twin;
    sf1_twin->twin = sf1;
    sf2->twin = sf2_twin;
    sf2_twin->twin = sf2;
}

// edge rolls to a face if singular; else rolls to a vertex 
void RollingMarkovModel::compute_edge_singularity_and_init_source_dir(){
    edge_is_singular = EdgeData<bool>(*mesh, false); // ~ is stable, by fwdSolver function names
    // initial assignment of real faces as SudoFaces
    root_sudo_face = HalfedgeData<SudoFace*>(*mesh);
    for (Edge e: mesh->edges()){
        Vertex next_vertex = forward_solver->next_rolling_vertex(e);
        
        if (next_vertex.getIndex() == INVALID_IND){ // singular edge
            edge_is_singular[e] = true;
            Halfedge he1 = e.halfedge(), 
                     he2 = e.halfedge().twin();
        }
        Vertex prev_vertex = e.otherVertex(next_vertex);
        Halfedge vec_field_aligned_he = (e.secondVertex() == next_vertex) ? e.halfedge() : e.halfedge().twin();
        // questionable; should I make the twin side null? or populate it?
        initiate_root_sudo_face(e);
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


