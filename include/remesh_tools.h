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


using namespace geometrycentral;
using namespace geometrycentral::surface;

void joint_remesh(ManifoldSurfaceMesh* mesh, VertexPositionGeometry *ref_geometry, VertexPositionGeometry *deformed_geometry, double target_edge_len);
bool split_only_remesh(ManifoldSurfaceMesh* mesh, VertexPositionGeometry *ref_geometry, VertexPositionGeometry *deformed_geometry, double target_edge_len);