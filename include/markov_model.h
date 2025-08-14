#pragma once

// #include "geometrycentral/surface/manifold_surface_mesh.h"
// #include "geometrycentral/surface/vertex_position_geometry.h"
// #include "geometrycentral/surface/surface_point.h"
#include "forward3D.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


// SudoFaces form a linked list, the root can be accessed from each edge 
// every SudoEdge formed by two SudoFaces is represented by the "first" SudoFace
class SudoFace{
    public:
        static size_t counter;
        const size_t index;
        Halfedge host_he;
        Vector3 normal;
        SudoFace *next_sudo_face, // on the host edge, just used for iteration 
                 *prev_sudo_face; // .. prev of first = first, and next of last = last
        // SudoFace *twin;
        SudoFace *source_sudo_face, // on neighboring edges
                 *sink_sudo_face;
        
        double EPS = 1e-8;
        
        // constructor
        SudoFace(Halfedge host_he_, Vector3 normal_, SudoFace *next_sudo_face_, SudoFace *prev_sudo_face_);
        // upon finding a new SudoFace
        SudoFace* split_sudo_edge(Vector3 new_normal);
};


class RollingMarkovModel {
    public:
        Forward3DSolver *forward_solver;

        // keeping G updated
        // Vector3 G;
        // assuming convexity
        // "same" pointer as the one in forward solver; here for easier access
        ManifoldSurfaceMesh* mesh;
        VertexPositionGeometry* geometry;
        
        double EPS = 1e-8;
        // constructors
        RollingMarkovModel(Forward3DSolver *forward_solver_);
        RollingMarkovModel(ManifoldSurfaceMesh *mesh_, VertexPositionGeometry* geometry_, Vector3 G_);

        // ---  one-time computable quantities ---

        // linked list data structure for smaller SudoArcs; needed for Markov chain edge surgery
        HalfedgeData<SudoFace*> root_sudo_face; // null on the sink side, potent on the source side (aligned with flow dir); if both not null, then we got a stabilizable edge (edge singularity)
        
        
        // start from sources and trace the vector field to split edge arcs and generate SudoEdges
        // initiate the recursive surgery
        std::list<Halfedge> termilar_hes; // bfs_list
        HalfedgeData<bool> he_processed;
        bool surgery_is_done = false;

        void initialize_pre_computes();

        // initialize sudo-face linked lists
        void initiate_root_sudo_face(Halfedge he);
        void init_root_sfs();
        
        // Split and surgery functions
        // handle single source HalfEdge/SudoEdge
        void process_halfedge(Halfedge he);
        void flow_he_to_he(Halfedge src, Halfedge dest);
        //includes a bunch of splits only on the dest SudoEdge
        void flow_sf_to_sf(SudoFace* src_sf1, SudoFace* dest_sf1);
        // building sudo faces and 
        void split_chain_edges_and_build_probability_pairs();
        
        // connected element pairs with non-zero probabilities 
        std::vector<std::pair<SudoFace*, SudoFace*>> sf_sf_pairs;
        std::vector<double> sf_sf_probs;
        std::vector<std::pair<SudoFace*, Face>> sf_face_pairs;
        std::vector<double> sf_face_probs;
        std::vector<std::pair<Vertex, SudoFace*>> vertex_sf_pairs;
        std::vector<double> vertex_sf_probs;
        // transition matrix for the whole Markov Chain
        SparseMatrix<double> transition_matrix;

        // clear pair/prob vectors
        void empty_prob_vectors();
        // sf to face pairs for non-trivial halfEdges
        void build_sf_face_pairs();
        
        // for probability pair debugging
        void print_prob_pairs();

        // build the transition matrix with Vertices, SudoEdges, Faces as nodes
        void build_transition_matrix();
        
        // check basic properties that should hold
        void check_transition_matrix();
        
        // probabilistic routes; only local routes (size-2 chains)
        double vertex_to_edge_prob(Vertex v, Edge e);
        double edge_to_edge_prob(Edge e1, Edge e2);
};