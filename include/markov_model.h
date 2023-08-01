#pragma once

// #include "geometrycentral/surface/manifold_surface_mesh.h"
// #include "geometrycentral/surface/vertex_position_geometry.h"
// #include "geometrycentral/surface/surface_point.h"
#include "forward3D.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

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

        // deterministic routes; assuming G is inside (positive mass), o.w. it won't be a DAG
        FaceData<Face> face_to_face; // might need to roll through a bunch of edges before getting to next face
        EdgeData<Face> edge_to_face;
        // probabilistic routes
        VertexData<std::vector<std::pair<Edge, double>>> vertex_to_edge_prob;
        EdgeData<std::vector<std::pair<Edge, double>>> edge_to_edge_prob;

};