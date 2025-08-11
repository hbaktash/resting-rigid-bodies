/************************************************************************
*
* ADOBE CONFIDENTIAL
* ___________________
*
* (Original header retained; consider replacing for publication.)
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
#include "polyscope/curve_network.h"

#include "forward3D.h"
#include "boundary_tools.h"
#include "geometry_utils.h"
#include "mesh_factory.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

static constexpr size_t ICOS_RES = 40;

// Holds only visualization state & tunable parameters (all static now)
class VisualUtils {
public:
  // polyscope handles
  inline static polyscope::SurfaceMesh* gm_sphere_mesh = nullptr;
  inline static polyscope::PointCloud* psG = nullptr;

  // radii
  inline static float stable_edge_radi = 0.007f;
  inline static float stablizable_edge_radi = 0.009f;
  inline static float both_edge_radi = 0.013f;
  inline static float gm_pt_radi = 0.005f;
  inline static float G_radi = 0.05f;
  inline static float arc_curve_radi = 0.002f;

  // colors
  inline static glm::vec3 stable_edge_color = glm::vec3({0.2f, 0.3f, 0.3f});
  inline static glm::vec3 stabilizable_edge_color = glm::vec3({0.05f, 0.05f, 0.05f});
  inline static glm::vec3 both_edge_color = glm::vec3({0.1f, 0.1f, 0.8f});
  inline static glm::vec3 patch_arc_fancy_color = glm::vec3({0.9f, 0.1f, 0.1f});

  // Gauss map config
  inline static double gm_distance = 2.0;
  inline static double gm_radi = 1.0;
  inline static Vector3 colored_shift = Vector3({gm_distance, gm_distance, 0.});
  inline static Vector3 center = Vector3({0., gm_distance, 0.});

  // flags
  inline static bool color_arcs = false;
  inline static bool gm_is_drawn = false;
  inline static bool draw_unstable_edge_arcs = true;
  inline static bool draw_stable_g_vec_for_unstable_edge_arcs = false;
  inline static bool show_hidden_stable_vertex_normals = false;

  // arc parameters
  inline static int arcs_seg_count = 15;
  inline static int arc_counter = 0;
};

// Free function interfaces (formerly member functions)
void draw_edge_arcs_on_gauss_map(Forward3DSolver* forwardSolver);
void draw_stable_vertices_on_gauss_map(Forward3DSolver* forwardSolver);
void draw_stable_face_normals_on_gauss_map(Forward3DSolver* forwardSolver);
void plot_height_function(Forward3DSolver* forwardSolver, ManifoldSurfaceMesh* sphere_mesh,
                          VertexPositionGeometry* sphere_geometry, bool plot_surface = true);
void draw_gauss_map(Forward3DSolver* forwardSolver, ManifoldSurfaceMesh* sphere_mesh,
                    VertexPositionGeometry* sphere_geometry);
void draw_guess_pc(Forward3DSolver* forwardSolver, std::vector<std::pair<size_t, size_t>> edge_inds,
                   std::vector<Vector3> boundary_normals);
void show_edge_equilibria_on_gauss_map(Forward3DSolver* forwardSolver);
void draw_G(Vector3 G);

void visualize_stable_vertices(Forward3DSolver* forwardSolver);
void visualize_edge_stability(Forward3DSolver* forwardSolver);
void visualize_face_stability(Forward3DSolver* forwardSolver);
void visualize_colored_polyhedra(Forward3DSolver* forwardSolver, FaceData<Vector3> face_colors);
void visualize_all_stable_orientations(Forward3DSolver* forwardSolver);

void draw_stable_patches_on_gauss_map(bool on_height_surface,
                                      BoundaryBuilder* bnd_builder,
                                      bool on_ambient_mesh);

void update_visuals(Forward3DSolver* tmp_solver, BoundaryBuilder* bnd_builder,
                    ManifoldSurfaceMesh* sphere_mesh, VertexPositionGeometry* sphere_geometry);

// Existing free utilities retained
void visualize_gauss_map(Forward3DSolver* forwardSolver);

// Arc helpers
void draw_arc_on_sphere(Vector3 p1, Vector3 p2, Vector3 center, double radius, size_t seg_count,
                        size_t edge_ind, double radi_scale = 1.,
                        glm::vec3 color = glm::vec3({-1.f, 0.f, 0.f}),
                        float arc_curve_radi = 0.003f);

void draw_arc_network_on_sphere(std::vector<std::pair<size_t, size_t>> edge_inds_,
                                std::vector<Vector3> positions_, Vector3 center, double radius,
                                size_t seg_count, std::string title, double radi_scale = 1.,
                                glm::vec3 color = glm::vec3({-1.f, 0.f, 0.f}),
                                float arc_curve_radi = 0.002f);

void draw_arc_network_on_lifted_suface(std::vector<std::pair<size_t, size_t>> edge_inds_,
                                       std::vector<Vector3> positions_,
                                       Forward3DSolver& forward_solver, Vector3 center,
                                       double radius, size_t seg_count, std::string title,
                                       glm::vec3 color = glm::vec3({-1.f, 0.f, 0.f}),
                                       float arc_curve_radi = 0.002f);

std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<Vector3>>
build_and_draw_stable_patches_on_gauss_map(BoundaryBuilder* boundary_builder,
                                           Vector3 center, double radius, size_t seg_count,
                                           bool on_height_surface = false,
                                           double arc_curve_radi = 0.004);

void draw_spherical_cone(std::vector<std::pair<size_t, size_t>> edges, std::vector<Vector3> poses,
                         Vector3 center, size_t seg_count, glm::vec3 color, std::string title);

std::tuple<ManifoldSurfaceMesh*, VertexPositionGeometry*>
make_cone_conforming_spherical_triangulation(BoundaryBuilder* boundary_builder,
    ManifoldSurfaceMesh* sphere_mesh, VertexPositionGeometry* sphere_geometry, Face f,
    Vector3 shift,
    std::vector<Vector3>& cone_poses, std::vector<std::pair<size_t, size_t>>& cone_edges,
    glm::vec3 cone_color,
    std::vector<Vector3>& prism_poses, std::vector<std::pair<size_t, size_t>>& prism_edges,
    glm::vec3 prism_color);

void visualize_face_solid_angle_vs_ms_complex(size_t f_ind, BoundaryBuilder* boundary_builder,
    ManifoldSurfaceMesh* sphere_mesh, VertexPositionGeometry* sphere_geometry);


void init_visuals(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry,
				  Forward3DSolver* forwardSolver, BoundaryBuilder* boundary_builder);


void draw_stable_patches_on_gauss_map(BoundaryBuilder* boundary_builder,
									  bool on_height_surface = false);

void draw_trail_on_gm(std::vector<Vector3> trail, glm::vec3 color, std::string name, double radi, bool color_gradient = false);

void visualize_quasi_static_drop_sequence(std::vector<Eigen::Matrix4d> transformation_matrices,
										 std::vector<Vector3> saved_snail_trail_refined,
										 Forward3DSolver* forwardSolver
										 );