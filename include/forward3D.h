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
    Face curr_face, next_face;
    Edge curr_edge, next_edge;
    Vertex curr_vertex, next_vertex;
    Vector3 curr_g_vec, next_g_vec;

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
    
    // 
    Edge vertex_to_edge(Vertex v, Vector3 curr_g_vec);
    
    void vertex_to_next(Vertex curr_v, Vector3 curr_g_vec);
    void edge_to_next(Edge curr_e, Vector3 curr_g_vec);
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