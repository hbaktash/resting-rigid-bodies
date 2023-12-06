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

// GC
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
// Polyscope
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"

#include "forward3D.h"
#include "boundary_tools.h"
#include "geometry_utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;



class VisualUtils{
    public:
        Forward3DSolver* forwardSolver;
        ManifoldSurfaceMesh* sphere_mesh;
        VertexPositionGeometry* sphere_geometry;
        polyscope::SurfaceMesh* gm_sphere_mesh;
        polyscope::PointCloud* psG;
        // polyscope::PointCloud *face_normals_p;
        float stable_edge_radi = 0.007,
              stablizable_edge_radi = 0.009,
              both_edge_radi = 0.013,
              gm_pt_radi = 0.01,
              G_radi = 0.1;

        
        glm::vec3 stable_edge_color = glm::vec3({0.2, 0.3, 0.3}),
                  stabilizable_edge_color = glm::vec3({0.05, 0.05, 0.05}),
                  // stabilizable_edge_color({0.1, 0.4, 0.6}),
                  both_edge_color = glm::vec3({0.1, 0.1, 0.8});


        // Gauss map stuff
        double gm_distance = 2.1,
               gm_radi = 1.0;
        Vector3 colored_shift = Vector3({gm_distance, gm_distance, 0.});
        // gm_shift = Vector3({0., gm_distance, 0.}),
                
        Vector3 center = Vector3({0., gm_distance, 0.});
        bool color_arcs = false, gm_is_drawn = false,
             draw_unstable_edge_arcs = true,
             draw_stable_g_vec_for_unstable_edge_arcs = false,
             show_hidden_stable_vertex_normals = true;

        // arc stuff
        // float arc_curve_radi = 0.0005;
        glm::vec3 patch_arc_fancy_color = glm::vec3({0.9,0.1,0.1});

        int arcs_seg_count = 12,
            arc_counter = 0;

        VisualUtils(){}
        VisualUtils(Forward3DSolver* forwardSolver){
                        this->forwardSolver = forwardSolver;
            // face_normal_vertex_gm_radi = 0.03;
        }

        void draw_edge_arcs_on_gauss_map();
        void draw_stable_vertices_on_gauss_map();
        void draw_stable_face_normals_on_gauss_map();
        void plot_height_function();
        void draw_gauss_map();
        void draw_guess_pc(std::vector<std::pair<size_t, size_t>> edge_inds, std::vector<Vector3> boundary_normals);
        void show_edge_equilibria_on_gauss_map();
        void draw_G();

        // polyhedra domain
        void visualize_stable_vertices();
        void visualize_edge_stability();
        void visualize_face_stability();
        void visualize_colored_polyhedra(FaceData<Vector3> face_colors);
        // void draw_arc_on_sphere();
        // void draw_arc_network_on_sphere();
        // void draw_arc_network_on_lifted_suface();
        // std::vector<Vector3> build_and_draw_stable_patches_on_gauss_map();
};


// arc stuff
// draw an arc connecting two points on the sphere; for Gauss map purposes
void draw_arc_on_sphere(Vector3 p1, Vector3 p2, Vector3 center, double radius, size_t seg_count, size_t edge_ind, polyscope::SurfaceMesh* hosting_psMesh, 
                        double radi_scale = 1., glm::vec3 color = glm::vec3({-1., 0, 0}), 
                        float arc_curve_radi = 0.003);


void draw_arc_network_on_sphere(std::vector<std::pair<size_t, size_t>> edge_inds_,
                                std::vector<Vector3> positions_,
                                Vector3 center, double radius, size_t seg_count, std::string title, 
                                polyscope::SurfaceMesh* hosting_psMesh, 
                                double radi_scale = 1., glm::vec3 color = glm::vec3({-1., 0, 0}), 
                                float arc_curve_radi = 0.003);

void draw_arc_network_on_lifted_suface(std::vector<std::pair<size_t, size_t>> edge_inds_,
                                       std::vector<Vector3> positions_,
                                       Forward3DSolver &forward_solver,
                                       Vector3 center, double radius, size_t seg_count, std::string title, 
                                       polyscope::SurfaceMesh* hosting_psMesh, 
                                       glm::vec3 color = glm::vec3({-1., 0, 0}), 
                                       float arc_curve_radi = 0.003);

std::pair<std::vector<std::pair<size_t, size_t>>,std::vector<Vector3>> 
    build_and_draw_stable_patches_on_gauss_map(BoundaryBuilder* boundary_builder, 
                                      polyscope::SurfaceMesh* hosting_psMesh,
                                      Vector3 center, double radius, size_t seg_count,
                                      bool on_height_surface = false);


