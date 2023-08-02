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
        Edge host_edge;
        Vector3 normal;
        SudoFace *next_sudo_face, // on the host edge, just used for iteration 
                 *prev_sudo_face; // .. prev of first = first, and next of last = last
        SudoFace *source_sudo_face, // on neighboring edges
                 *sink_sudo_face;

        // constructor
        SudoFace(Edge host_edge_, Vector3 normal_, SudoFace *next_sudo_face_, SudoFace *prev_sudo_face_);
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

        // for building the actual Markov model
        VertexData<bool> vertex_stabilizablity;
        HalfedgeData<SudoFace*> first_sudo_face; // null on the sink side, potent on the source side; if both not null, then we got a stabilizable edge (edge singularity)

        // deterministic routes; assuming G is inside (positive mass), o.w. it won't be a DAG (will have loops)
        FaceData<Face> face_to_face; // might need to roll through a bunch of edges before getting to next face
        EdgeData<Face> edge_to_face;
        // probabilistic routes
        VertexData<std::vector<std::pair<Edge, double>>> vertex_to_edge_prob;
        EdgeData<std::vector<std::pair<Edge, double>>> edge_to_edge_prob;

};