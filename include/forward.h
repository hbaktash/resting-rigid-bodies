#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class ForwardSolver {

  public:
    VertexData<Vertex> next;
    VertexData<Vertex> prev;
    //2D geometry; z = 0 always
    Vector3 G;
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;

    // state: just using a pair
    std::pair<Vertex, Vertex> curr_state;
    std::string current_state = "vertex"; //w.p. 1
    // constructors
    ForwardSolver() {}
    ForwardSolver(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo);

    void center_geometry();
    void align_geometry(Vector2 dir);
    // for one instance
    void find_contact(Vector2 initial_ori);
    void next_state();

    void integrate(double h);
};