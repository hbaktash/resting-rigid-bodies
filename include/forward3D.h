#pragma once

#include "arc_algebra.h"
#include "geometry_utils.h"
#include "convex_hull.h"

// #include "geometrycentral/surface/manifold_surface_mesh.h"
// #include "geometrycentral/surface/vertex_position_geometry.h"
// #include "geometrycentral/surface/surface_point.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


class Forward3DSolver {
  private:
    Vector3 G;
  public:

    double volume = -1.;
    // input goemetry
    ManifoldSurfaceMesh* inputMesh;
    VertexPositionGeometry* inputGeometry;
    // convex hull of input goemetry (a polyhedra) and center of mass
    ManifoldSurfaceMesh* hullMesh;
    VertexPositionGeometry* hullGeometry;

    // hull <-> mesh mapping
    VertexData<size_t> org_hull_indices; //hull vertex indices in the original mesh
    VertexData<size_t> on_hull_index; // index of an interior vertex on the hull; INVALID_IND if not on hull
    Vector<size_t> hull_indices, interior_indices;
    void trivial_initialize_index_trackers();

    // current state; for simulation 
    Vertex curr_v;
    Edge curr_e;
    Face curr_f;
    Vector3 curr_g_vec;
    bool stable_state = false;
    // for later vector field display
    Vector3 initial_roll_dir;

    // for debugging purposes
    Vector3 tmp_test_vec;

    // constructors
    Forward3DSolver() {}
    // TODO; make constructor for non_convex input
    Forward3DSolver(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, Vector3 G, bool concave_input = true);
    Forward3DSolver(Eigen::MatrixX3d point_cloud, Eigen::Vector3d G, bool is_convex);
    // setter/getter for G
    void set_G(Vector3 new_G);
    void set_uniform_G();
    Vector3 get_G();

    // Hull related updates
    // bool first_hull = true;
    // void update_hull_index_arrays();
    // void update_convex_hull(bool with_projection = false);
    // void update_hull_points_correspondence(VertexData<Vector3> new_hull_points, VertexData<Vector3> old_points);

    // initialize state
    void initialize_state(Vertex curr_v, Edge curr_e, Face curr_f, Vector3 curr_g_vec);
    
    // heigh for a given normal that depends on contact point
    double height_function(Vector3 ground_normal);
    
    // find initial contact vertex
    void find_contact(Vector3 initial_ori);

    // build next element DP stuff
    void compute_vertex_probabilities();
    // stabilizable := the normal from G can touch the ground
    // stable := the normal from G falls withing the element
    bool vertex_is_stablizable(Vertex v);
    Vertex next_rolling_vertex(Edge e); // INVALID index if edge is singular; i.e. stable
    bool edge_is_stable(Edge e);
    bool edge_is_stablizable(Edge e);
    bool face_is_stable(Face f);
  
    
    // pre-compute containers
    // whether the stable point can be reached or not
    VertexData<bool> vertex_is_stabilizable;
    // stable normal of a vertex; even if not reachable
    VertexData<Vector3> vertex_stable_normal;
    VertexData<double> vertex_gaussian_curvature;
    VertexData<double> vertex_probabilities;
    // edge rolls to a face if singular; else rolls to a vertex
    EdgeData<bool> edge_is_singular;
    EdgeData<Vertex> edge_next_vertex;
    // is 0 , if the normal is unreachable, or doesnt fall on the edge (edge too short)
    EdgeData<Vector3> edge_stable_normal;
    // deterministic routes; face to next face
    FaceData<Face> face_next_face; // might need to roll through a bunch of edges before getting to next face
    FaceData<Face> face_last_face; // keep rolling till u hit a stable face
        
    // pre-compute functions
    void compute_vertex_stabilizablity();
    void compute_edge_stable_normals();
    void compute_vertex_gaussian_curvatures();
    void build_face_next_faces();
    // not basic; calls prev function
    void build_face_last_faces();

    // just calls basic precomputes
    void initialize_pre_computes();    
    void print_precomputes();
    // Rigid Simulation stuff
    Edge vertex_to_edge(Vertex v);
    void vertex_to_next(Vertex v);
    void edge_to_next(Edge e);
    void face_to_next(Face f);
    
    // rigid simulation
    std::vector<Vector3> translation_log;
    std::vector<Vertex> vertex_log;
    std::vector<Edge> edge_log;
    std::vector<Face> face_log;
    void next_state(bool verbose = false);
    
    std::vector<Vector3> snail_trail_log(Vector3 initial_orientation);
    // // empirical sampling
    Face final_touching_face(Vector3 initial_ori);
    // void empirically_build_probabilities(int sample_count);

};