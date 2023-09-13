/************************************************************************
*
* ADOBE CONFIDENTIAL
* ___________________
*
* Copyright [first year code created] Adobe
* All Rights Reserved.
*
* NOTICE: All information contained herein is, and remains
* the property of Adobe and its suppliers, if any. The intellectual
* and technical concepts contained herein are proprietary to Adobe
* and its suppliers and are protected by all applicable intellectual
* property laws, including trade secret and copyright laws.
* Dissemination of this information or reproduction of this material
* is strictly forbidden unless prior written permission is obtained
* from Adobe.
*************************************************************************
*/
#pragma once

#include "arc_algebra.h"
// #include "geometrycentral/surface/manifold_surface_mesh.h"
// #include "geometrycentral/surface/vertex_position_geometry.h"
// #include "geometrycentral/surface/surface_point.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;

Vector3 project_on_plane(Vector3 p, Vector3 offset, Vector3 normal);
Vector3 point_to_segment_normal(Vector3 P, Vector3 A, Vector3 B);

class Forward3DSolver {

  public:
    //input convex goemetry (a polyhedra) and center of mass
    Vector3 G;
    ManifoldSurfaceMesh* hullMesh;
    VertexPositionGeometry* hullGeometry;

    // current state; for simulation 
    Vertex curr_v;
    Edge curr_e;
    Face curr_f;
    Vector3 curr_g_vec;
    bool stable_state = false;
    // for later vector field display
    Vector3 initial_roll_dir;

    // just gaussian curvature of a vertex
    VertexData<double> vertex_probabilities;

    // for debugging purposes
    Vector3 tmp_test_vec;

    // constructors
    Forward3DSolver() {}
    Forward3DSolver(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo, Vector3 G);

    // initialize state
    void initialize_state(Vertex curr_v, Edge curr_e, Face curr_f, Vector3 curr_g_vec);
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
    // edge rolls to a face if singular; else rolls to a vertex
    EdgeData<bool> edge_is_singular;
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