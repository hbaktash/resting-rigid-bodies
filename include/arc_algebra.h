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

// #include <autodiff/reverse/var/eigen.hpp>
// #include <autodiff/reverse/var.hpp>

using namespace geometrycentral;
using namespace geometrycentral::surface;


bool is_on_arc_segment(Vector3 nomral, Vector3 A, Vector3 B);

// A-portion
double arc_portion(Vector3 mid_point, Vector3 A, Vector3 B);

double patch_area(Vector3 A, Vector3 B, Vector3 C, Vector3 D);

// 0 if not within AB, else intersection normal with norm 1
Vector3 intersect_arc_ray_with_arc(Vector3 R1, Vector3 R2, Vector3 A, Vector3 B, bool &sign_change);

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


// not so arc
Vector3 project_on_plane(Vector3 p, Vector3 offset, Vector3 normal);

Vector3 point_to_segment_normal(Vector3 P, Vector3 A, Vector3 B);


// autodiff stuff; should use templates at some point
// autodiff::Vector3var point_to_segment_normal_ad(autodiff::MatrixX3var &poses, autodiff::Vector3var &G, Edge e);
// autodiff::Vector3var face_normal_ad(autodiff::MatrixX3var &poses, Face f);
// autodiff::Vector3var intersect_arc_ray_with_arc_ad(autodiff::MatrixX3var &poses, autodiff::Vector3var &G, Vertex v,
//                                                    autodiff::Vector3var &R2, autodiff::Vector3var &A, autodiff::Vector3var &B,
//                                                    bool sign_change);
// autodiff::var triangle_patch_area_on_sphere_ad(autodiff::Vector3var &A, autodiff::Vector3var &B, autodiff::Vector3var &C);
