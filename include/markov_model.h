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
        Halfedge host_he;
        Vector3 normal;
        SudoFace *next_sudo_face, // on the host edge, just used for iteration 
                 *prev_sudo_face; // .. prev of first = first, and next of last = last
        SudoFace *twin;
        SudoFace *source_sudo_face, // on neighboring edges
                 *sink_sudo_face;

        // constructor
        SudoFace(Halfedge host_he_, Vector3 normal_, SudoFace *next_sudo_face_, SudoFace *prev_sudo_face_);
        // upon finding a new SudoFace
        SudoFace* split_sudo_edge(Vector3 new_normal);
};


class RollingMarkovModel {
    public:
        Forward3DSolver *forward_solver;

        // keeping G updated
        Vector3 G;
        // assuming convexity
        // "same" pointer as the one in forward solver; here for easier access
        ManifoldSurfaceMesh* mesh;
        VertexPositionGeometry* geometry;

        // constructors
        RollingMarkovModel(Forward3DSolver *forward_solver_);
        RollingMarkovModel(ManifoldSurfaceMesh* mesh_, VertexPositionGeometry* geometry_, Vector3 G_);

        // ---  one-time computable quantities ---

        // whether the stable point can be reached or not
        VertexData<bool> vertex_stabilizablity; 
        VertexData<Vector3> vertex_stable_normal; 
        void compute_vertex_stabilizablity();
        // either reachable or un-reachable; inferable from _vertex_stabilizability_
        // edge rolls to a face if singular; else rolls to a vertex
        EdgeData<bool> edge_is_singular; 
        void compute_edge_singularity_and_init_source_dir();
        // is 0 , if the normal is unreachable, or doesnt fall on the edge (edge too short)
        EdgeData<Vector3> edge_stable_normal; 
        void compute_edge_stable_normal();


        // deviding arcs t SudoEdges; Markov chain edge surgery
        HalfedgeData<SudoFace*> root_sudo_face; // trivial (or null??) on the sink side, potent on the source side (aligned with flow dir); if both not null, then we got a stabilizable edge (edge singularity)
        void initiate_root_sudo_face(Edge e); // deals with both he's

        // deterministic routes; assuming G is inside (positive mass), o.w. it won't be a DAG (will have loops)
        FaceData<Face> face_to_face; // might need to roll through a bunch of edges before getting to next face
        EdgeData<Face> edge_to_face;
        // probabilistic routes; only local routes (size-2 chains)
        double vertex_to_edge_prob(Vertex v, Edge e);
        double edge_to_edge_prob(Edge e1, Edge e2);

};