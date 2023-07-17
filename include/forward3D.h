#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_point.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;


Vector3 project_on_plane(Vector3 p, Vector3 offset, Vector3 normal);

class Forward3DSolver {

  public:
    //input goemetry (a polyhedra) and center of mass
    Vector3 G;
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;

    // convex hull geometry; another polyhedra
    ManifoldSurfaceMesh* hullMesh;
    VertexPositionGeometry* hullGeometry;

    // current state; if actually trying to simulate
    Vertex curr_v1, curr_v2, curr_v3;
    Vector3 curr_g_vec;

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
    Edge vertex_to_edge(Vertex v);
    
    void vertex_to_next(Vertex curr_v);
    void edge_to_next(Edge curr_e);
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