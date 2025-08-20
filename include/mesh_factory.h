#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/remeshing.h"

#include "geometry_utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

Vector3 spherical_to_xyz(double r, double phi, double theta);
Vector3 cylindrical_to_xyz(double h, double r, double theta);

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> 
generate_polyhedra(std::string poly_str);

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> 
generate_11_sided_polyhedron(std::string type);

void preprocess_mesh(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry,
                     bool triangulate = false, bool do_remesh = false, double remesh_edge_scale = 3., bool normalize = true);

// sphere tesselation generation
std::vector<Vector3> generate_normals_icosahedral(int resolution);

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> 
generate_pointy_prism(size_t n);

std::tuple<ManifoldSurfaceMesh*, VertexPositionGeometry*>
icosahedron();