#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


double signed_volume(Vector3 A, Vector3 B, Vector3 C, Vector3 D);


std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> 
convex_hull(ManifoldSurfaceMesh &input_mesh, VertexPositionGeometry &input_geometry);


// assumes uniform mass density; TODO how to handle non-uniform?
Vector3 find_center_of_mass(ManifoldSurfaceMesh &input_mesh, VertexPositionGeometry &input_geometry);

void center_and_normalize(ManifoldSurfaceMesh* input_mesh, VertexPositionGeometry* input_geometry);