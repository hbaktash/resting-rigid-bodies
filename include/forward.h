#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class ForwardSolver {

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
    EdgeData<Edge> next_falling_edge;
    VertexData<double> vertex_probabilities;
    EdgeData<double> initial_edge_probabilities;
    EdgeData<double> final_edge_probabilities;


    // constructors
    ForwardSolver() {}
    ForwardSolver(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo,
                  Vector3 G);

    // find initial contact vertex
    void find_contact(Vector3 initial_ori);

    // build next edge DP stuff
    void build_next_edge_tracer();
    void compute_vertex_probabilities();
    void compute_initial_edge_probabilities();
    void compute_final_edge_probabilities();

    // lazy polygon navigation
    Edge other_edge(Edge curr_e, Vertex tip_v);

    // maybe later for fancier visualization
    void center_geometry();
    void align_geometry(Vector3 dir);
    void next_state();
    

};