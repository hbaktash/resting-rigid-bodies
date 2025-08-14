#pragma once

#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/Qhull.h"
#include "libqhull/qhull_a.h"

#include "utils.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using orgQhull::Qhull;
using orgQhull::QhullFacet;
using orgQhull::QhullVertex;
using namespace geometrycentral;
using namespace geometrycentral::surface;


// std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> 
std::tuple<std::vector<std::vector<size_t>>, std::vector<size_t>, std::vector<Vector3>>
get_convex_hull(Eigen::MatrixX3d point_cloud);

std::tuple<ManifoldSurfaceMesh*, VertexPositionGeometry*>
get_convex_hull_mesh(Eigen::MatrixX3d point_cloud);

std::tuple<std::vector<std::vector<size_t>>, std::vector<size_t>, std::vector<Vector3>>
get_convex_hull(std::vector<Vector3> point_set);

std::tuple<std::vector<std::vector<size_t>>, std::vector<size_t>, std::vector<Vector3>>
get_convex_hull(VertexData<Vector3> point_set);


std::tuple<ManifoldSurfaceMesh*, VertexPositionGeometry*>
get_convex_hull_mesh(std::vector<Vector3> point_set);

std::tuple<ManifoldSurfaceMesh*, VertexPositionGeometry*>
get_convex_hull_mesh(VertexData<Vector3> point_set);

std::tuple<ManifoldSurfaceMesh*, VertexPositionGeometry*> 
get_mesh_for_convex_set(Eigen::MatrixX3d convex_point_cloud);

Vector3 project_back_into_hull(VertexPositionGeometry *hull_geometry, Vector3 p);