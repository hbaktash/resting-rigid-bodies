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
                     bool triangulate = false, bool do_remesh = false, double remesh_edge_scale = 3.);

// sphere tesselation generation
std::vector<Vector3> generate_normals_icosahedral(int resolution);