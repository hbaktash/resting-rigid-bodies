#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/remeshing.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;

void joint_remesh(ManifoldSurfaceMesh* mesh, VertexPositionGeometry *ref_geometry, VertexPositionGeometry *deformed_geometry, double target_edge_len);
bool split_only_remesh(ManifoldSurfaceMesh* mesh, VertexPositionGeometry *ref_geometry, VertexPositionGeometry *deformed_geometry, double target_edge_len);