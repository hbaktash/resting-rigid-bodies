#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_point.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;


class StablePoint {
  public:
    Vector3 position;
    Vertex host_vertex = Vertex();
    Edge host_edge = Edge();
    std::vector<StablePoint> adj_stable_points;

    // constructors
    StablePoint(){}
    StablePoint(Vertex host_v, Vector3 position);
    StablePoint(Edge host_e, Vector3 position);

    // ipp stuff
    bool is_vertex();
    bool is_edge();

};


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

    // current state; simulation 
    Vertex curr_v;
    Edge curr_e;
    Face curr_f;
    Vector3 curr_g_vec;
    bool stable_state = false;
    // for later display
    Vector3 initial_roll_dir;

    // DP stuff; implicitly using a DAG
    FaceData<Face> next_falling_face;
    VertexData<double> vertex_probabilities;
    
    // for debugging purposes
    Vector3 tmp_test_vec;

    // constructors
    Forward3DSolver() {}
    Forward3DSolver(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, Vector3 G);

    // find initial contact vertex
    void find_contact(Vector3 initial_ori);

    // build next element DP stuff
    void compute_vertex_probabilities();
    // stabilizable := the normal from G can touch the ground
    // stable := the normal from G falls withing the element
    bool vertex_is_stablizable(Vertex v);
    bool edge_is_stable(Edge e);
    bool edge_is_stablizable(Edge e);
    bool face_is_stable(Face f);
  
    // Rigid Simulation stuff
    Edge vertex_to_edge(Vertex v);
    void vertex_to_next(Vertex v);
    void edge_to_next(Edge e);
    void face_to_next(Face f);
    
    // rigid simulation
    void next_state();
    
    // // empirical sampling
    Face final_touching_face(Vector3 initial_ori);
    // void empirically_build_probabilities(int sample_count);

};