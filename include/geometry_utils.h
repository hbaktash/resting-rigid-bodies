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

using namespace geometrycentral;
using namespace geometrycentral::surface;


double signed_volume(Vector3 A, Vector3 B, Vector3 C, Vector3 D);


std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> 
convex_hull(ManifoldSurfaceMesh &input_mesh, VertexPositionGeometry &input_geometry);


bool G_is_inside(ManifoldSurfaceMesh &mesh, VertexPositionGeometry &geometry, Vector3 G);

// Gaussian curvature for a general polygonal case
bool gaussian_curvature(Vertex v, VertexPositionGeometry &geometry);


// assumes uniform mass density; TODO how to handle non-uniform?
Vector3 find_center_of_mass(ManifoldSurfaceMesh &input_mesh, VertexPositionGeometry &input_geometry);

void center_and_normalize(ManifoldSurfaceMesh* input_mesh, VertexPositionGeometry* input_geometry);