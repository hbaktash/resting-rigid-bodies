#pragma once

// GC
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
// Polyscope
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


// draw an arc connecting two points on the sphere; for Gauss map purposes
void draw_arc_on_sphere(Vector3 p1, Vector3 p2, Vector3 center, double radius, size_t seg_count, size_t edge_ind, polyscope::SurfaceMesh* hosting_psMesh, 
                        double radi_scale = 1., glm::vec3 color = glm::vec3({-1., 0, 0}), 
                        float arc_curve_radi = 0.01, 
                        glm::vec3 default_arc_color = glm::vec3({0.05, 0.05, 0.05}),
                        glm::vec3 default_patch_arc_color = glm::vec3({0.05, 0.5, 0.5}),
                        glm::vec3 snail_trail_color = glm::vec3({0.8,0.8,0.2}));