#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class Forward3DSolver {

  public:
    //2D geometry; z = 0 always

    //input goemetry (a polygon) and center of mass
    Vector3 G;
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;

    // convex hull geometry; another polygon
    ManifoldSurfaceMesh* hullMesh;
    VertexPositionGeometry* hullGeometry;

    // current state; if actually trying to simulate
    std::pair<Vertex, Vertex> curr_state;
    Vector3 current_g_vec;

    // DP stuff; implicitly using a DAG
    FaceData<Face> next_falling_face;
    VertexData<double> vertex_probabilities;
    
    // for debugging purposes
    Vector3 tmp_test_vec;

    // constructors
    Forward3DSolver() {}
    Forward3DSolver(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, Vector3 G);

    // find initial contact vertex
    // void find_contact(Vector3 initial_ori);

    // build next element DP stuff
    void compute_vertex_probabilities();
    // stabilizable := the normal from G can touch the ground
    // stable := the normal from G falls withing the element
    bool vertex_is_stablizable(Vertex v);
    bool edge_is_stable(Edge e);
    bool edge_is_stablizable(Edge e);
    bool face_is_stable(Face f);
    // void compute_initial_edge_probabilities();
    
    // void compute_final_edge_probabilities();

    // // lazy polygon navigation
    // Edge other_edge(Edge curr_e, Vertex tip_v);

    // // possible fancy visualization
    // void center_geometry();
    // void align_geometry(Vector3 dir);
    // // rigid simulation
    // void next_state();
    
    // // empirical sampling
    // Edge simulate_toss(Vector3 initial_ori);
    // void empirically_build_probabilities(int sample_count);


};