#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

Vector3 spherical_to_xyz(double r, double phi, double theta);
Vector3 cylindrical_to_xyz(double h, double r, double theta);

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> 
generate_polyhedra(std::string poly_str);