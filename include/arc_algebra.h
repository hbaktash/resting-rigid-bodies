#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;

// 0 if not within AB, else intersection normal with norm 1
Vector3 intersect_arc_ray_with_arc(Vector3 R1, Vector3 R2, Vector3 A, Vector3 B);

// 
// Vector3 intersect_arc_ray_with_arc(Vector3 R1, Vector3 R2, Vector3 A, Vector3 B);