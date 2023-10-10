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


bool is_on_arc_segment(Vector3 nomral, Vector3 A, Vector3 B);

// A-portion
double arc_portion(Vector3 mid_point, Vector3 A, Vector3 B);

double patch_area(Vector3 A, Vector3 B, Vector3 C, Vector3 D);

// 0 if not within AB, else intersection normal with norm 1
Vector3 intersect_arc_ray_with_arc(Vector3 R1, Vector3 R2, Vector3 A, Vector3 B);

// 
// Vector3 intersect_arc_ray_with_arc(Vector3 R1, Vector3 R2, Vector3 A, Vector3 B);


// Gaussian curvature for a general polygonal case
double gaussian_curvature(Vertex v, VertexPositionGeometry &geometry);

// area of a triangular patch on sphere
double triangle_patch_area_on_sphere(Vector3 A, Vector3 B, Vector3 C);

// angle between two spherical arcs; <P1AP2
double angle_on_sphere(Vector3 P1, Vector3 A, Vector3 P2);

// angle between two spherical arcs; <P1AP2
double angle_on_sphere(Vector3 P1, Vector3 A, Vector3 P2);

