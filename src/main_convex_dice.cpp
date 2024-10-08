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
#include <stan/math.hpp>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "args.hxx"
#include "imgui.h"

#include "coloring.h"
#include "forward3D.h"
#include "mesh_factory.h"
// #include "geometry_utils.h" //; in fwd solver 
#include "visual_utils.h"
#include "inv_design.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


#define ANSI_FG_MAGENTA "\x1b[35m"
#define ANSI_FG_YELLOW "\x1b[33m"
#define ANSI_FG_GREEN "\x1b[32m"
#define ANSI_FG_WHITE "\x1b[37m"
#define ANSI_FG_RED "\x1b[31m"
#define ANSI_RESET "\x1b[0m"


// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr, cv_mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr, cv_geometry_ptr;
ManifoldSurfaceMesh *mesh, *convex_to_fill_mesh;
VertexPositionGeometry *geometry, *convex_to_fill_geometry, *deformed_geometry;
bool load_hull_from_file = true;
Vector3 G, // center of Mass
        pre_deform_G({-10.,-10.,-10.}),
        post_deform_G({-10.,-10.,-10.}),
        initial_g_vec({0,-1,0}),
        default_face_color({0.99,0.99,0.99}),
        curr_face_color({0.1,0.87,0.1});

Forward3DSolver* forwardSolver;
  

VisualUtils vis_utils;

ManifoldSurfaceMesh* sphere_mesh;
VertexPositionGeometry* sphere_geometry;
        

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psInputMesh, *psHullMesh;

polyscope::PointCloud *raster_pc;

float pt_cloud_radi_scale = 0.1,
      G_r = 0.,
      G_theta = 0.,
      G_phi = 0.;


bool gm_is_drawn;
bool color_arcs = true,
     draw_unstable_edge_arcs = true,
     show_hidden_stable_vertex_normals = false;

// raster image stuff
int sample_count = 8000; 
FaceData<Vector3> face_colors;
bool recolor_faces = true,
     real_time_raster = false,
     draw_point_cloud = true,
     normalize_vecF = true;
glm::vec3 vecF_color({0.1, 0.1, 0.1});

// boundary stuff
BoundaryBuilder *boundary_builder;
polyscope::PointCloud *boundary_normals_pc;
polyscope::SurfaceMesh *dummy_psMesh_for_regions, *dummy_psMesh_for_height_surface;
bool draw_boundary_patches = false; // TODO
bool test_guess = false;
bool compute_global_G_effect = true,
     deform_after = true,
     frozen_G = true,
     structured_opt = true,
     dynamic_remesh = true,
     use_reg = false,
     use_QP_solver = true,
     always_update_structure = true,
     with_hull_projection = false,
     curvature_weighted_CP = false;
bool do_remesh = false;
float scale_for_remesh = 1.003;
polyscope::PointCloud *test_pc;

int fair_sides_count = 6; // for optimization
bool do_sobolev_dice_grads = true,
     use_autodiff_for_dice_grad = true;
float sobolev_lambda = 2.,
      sobolev_lambda_decay = 0.95,
      dice_energy_step = 0.01;
int sobolev_p = 2;
// optimization stuff

float dice_search_decay = 0.95;


// example choice
std::vector<std::string> all_input_names = {std::string("triangular"), std::string("circus tent")}; //{std::string("tet"), std::string("tet2"), std::string("cube"), std::string("tilted cube"), std::string("dodecahedron"), std::string("Conway spiral 4"), std::string("oloid")};
std::string input_name = "triangular";


void draw_stable_patches_on_gauss_map(bool on_height_surface = false, 
                                      BoundaryBuilder *bnd_builder = boundary_builder,
                                      Forward3DSolver *tmp_solver = forwardSolver){
  auto net_pair = build_and_draw_stable_patches_on_gauss_map(bnd_builder, vis_utils.center, vis_utils.gm_radi, vis_utils.arcs_seg_count, on_height_surface);
  if (test_guess)
    vis_utils.draw_guess_pc(net_pair.first, net_pair.second);
}


void visualize_gauss_map(){
  std::unique_ptr<ManifoldSurfaceMesh> sphere_mesh_ptr;
  std::unique_ptr<VertexPositionGeometry> sphere_geometry_ptr;
  std::tie(sphere_mesh_ptr, sphere_geometry_ptr) = generate_polyhedra("sphere");
  sphere_mesh = sphere_mesh_ptr.release();
  sphere_geometry = sphere_geometry_ptr.release();
  
  //update vis utils
  vis_utils.sphere_geometry = sphere_geometry;
  vis_utils.sphere_mesh = sphere_mesh;

  vis_utils.draw_gauss_map();
  gm_is_drawn = true;
}


void init_visuals(){
  // Register the mesh with polyscope
  psInputMesh = polyscope::registerSurfaceMesh(
    "init input mesh",
    geometry->inputVertexPositions, mesh->getFaceVertexList(),
    polyscopePermutations(*mesh));
  psInputMesh->setTransparency(0.75);
  psInputMesh->setEnabled(true);
  psHullMesh = polyscope::registerSurfaceMesh(
    "init hull mesh",
    forwardSolver->hullGeometry->inputVertexPositions, forwardSolver->hullMesh->getFaceVertexList(),
    polyscopePermutations(*forwardSolver->hullMesh));
  psHullMesh->setEnabled(false);
  vis_utils.forwardSolver = forwardSolver;
  vis_utils.draw_G();
}

void update_solver(){
  //assuming convex input here
  forwardSolver = new Forward3DSolver(mesh, geometry, G);
  forwardSolver->set_uniform_G();
  forwardSolver->initialize_pre_computes();
  boundary_builder = new BoundaryBuilder(forwardSolver);
  boundary_builder->build_boundary_normals();
  vis_utils.forwardSolver = forwardSolver;
}


void generate_polyhedron_example(std::string poly_str, bool triangulate = false){
    // readManifoldSurfaceMesh()
//   std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(poly_str);
  std::tie(mesh_ptr, geometry_ptr) = generate_11_sided_polyhedron(poly_str);
  mesh = mesh_ptr.release();
  geometry = geometry_ptr.release();
  preprocess_mesh(mesh, geometry, triangulate, false, scale_for_remesh);
}



void update_visuals(Forward3DSolver *tmp_solver = nullptr, BoundaryBuilder *bnd_builder = boundary_builder){
  vis_utils.forwardSolver = tmp_solver;
  vis_utils.draw_gauss_map();
  vis_utils.draw_G();
  vis_utils.plot_height_function();
  draw_stable_patches_on_gauss_map(false, bnd_builder);
}

void initialize_state(std::string input_name){
    generate_polyhedron_example(input_name);
    update_solver();
    init_visuals();
    visualize_gauss_map();//
    boundary_builder->print_area_of_boundary_loops();
    update_visuals(forwardSolver, boundary_builder);
}

void get_dice_energy_grads(Eigen::MatrixX3d hull_positions, Eigen::Vector3d G_vec,
                           Eigen::MatrixX3d &df_dv, Eigen::Vector3d &df_dG, double &dice_energy,
                           bool use_autodiff, std::string policy, int fair_sides){
    if (!use_autodiff){
      Forward3DSolver fwd_solver(hull_positions, G_vec);
      fwd_solver.initialize_pre_computes();
      printf("   initialized for diffs\n");
      BoundaryBuilder *bnd_builder = new BoundaryBuilder(&fwd_solver);
      bnd_builder->build_boundary_normals();
      InverseSolver *inv_solver = new InverseSolver(bnd_builder);
      inv_solver->find_uni_mass_d_pf_dv(false, frozen_G);
      // distribution is set internally here
      VertexData<Vector3> approx_dice_energy_grads = inv_solver->find_uni_mass_total_vertex_grads(fair_sides_count,
                                                                                                  0.01);
      inv_solver->find_d_pf_d_Gs(false);
      Vector3 G_grad = inv_solver->find_total_g_grad();
      // convert to eigen
      dice_energy = bnd_builder->get_fair_dice_energy(fair_sides);
      df_dv = vertex_data_to_matrix(approx_dice_energy_grads);
      df_dG = vec32vec(G_grad);
    }
    else {
      Forward3DSolver fwd_solver(hull_positions, G_vec);
      fwd_solver.initialize_pre_computes();
      Eigen::MatrixX3d hull_positions = vertex_data_to_matrix(fwd_solver.hullGeometry->inputVertexPositions); 
      auto dice_energy_lambda = [&] <typename Scalar> (const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &hull_poses_G_append_vec) -> Scalar {
        // decompose flat vector to positions and center of mass; G is the last 3 elements
        Eigen::Vector3<Scalar> G_eigen = hull_poses_G_append_vec.tail(3);
        size_t flat_n = hull_poses_G_append_vec.rows();
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 3> > hull_poses(hull_poses_G_append_vec.head(flat_n-3).data(), flat_n/3 - 1, 3);
        return BoundaryBuilder::dice_energy<Scalar>(hull_poses, G_eigen, fwd_solver, fair_sides);
      };
      Eigen::VectorXd hull_poses_vec = hull_positions.reshaped();
      Eigen::VectorXd hull_poses_and_G_vec(hull_poses_vec.size() + 3);
      hull_poses_and_G_vec << hull_poses_vec, G_vec;
      
      Eigen::VectorXd dfdU_vec;
      double dice_e;
      stan::math::gradient(dice_energy_lambda, hull_poses_and_G_vec, dice_e, dfdU_vec);
      dice_energy = dice_e;
      size_t flat_n = dfdU_vec.rows();
      df_dG = dfdU_vec.tail(3);
      Eigen::Map<Eigen::MatrixXd> dfdV(dfdU_vec.head(flat_n-3).data(), flat_n/3 - 1, 3);
      df_dv = dfdV;
    }
    
}

void dice_energy_opt(size_t step_count = 1){
  Forward3DSolver tmp_solver(mesh, geometry, G, true);
  tmp_solver.set_uniform_G();
  G = tmp_solver.get_G();
  tmp_solver.initialize_pre_computes();
  Eigen::Vector3d G_vec{G.x, G.y, G.z};
  Eigen::MatrixX3d hull_positions = vertex_data_to_matrix(tmp_solver.hullGeometry->inputVertexPositions);
  std::cout << "initial G: " << G_vec << "\n";
  
  // BoundaryBuilder tmp_bnd_builder(&tmp_solver);
  
  double current_sobolev_lambda = sobolev_lambda;

  for (size_t iter = 0; iter < step_count; iter++){
    double dice_e;
    Eigen::Vector3d dfdG;
    Eigen::MatrixX3d dfdV;
    
    printf("getting grads\n");
    get_dice_energy_grads(hull_positions, G_vec, dfdV, dfdG, dice_e, use_autodiff_for_dice_grad, 
                          "[policy]", fair_sides_count);
    // update_visuals_with_G(&tmp_solver, &tmp_bnd_builder);
    std::cout << ANSI_FG_YELLOW << " DE at iter "<< iter << " f: " << dice_e << ANSI_RESET << "\n";
    auto curr_hull_psmesh = polyscope::registerSurfaceMesh("current hull", tmp_solver.hullGeometry->inputVertexPositions, tmp_solver.hullMesh->getFaceVertexList());
    curr_hull_psmesh->setSurfaceColor({0.1,0.9,0.1})->setEdgeWidth(2.)->setTransparency(0.7)->setEnabled(true);
    // show grads
    curr_hull_psmesh->addVertexVectorQuantity("ad total grads", -1.*dfdV)->setEnabled(true);
    // diffused grads
    if (do_sobolev_dice_grads){    
      current_sobolev_lambda *= sobolev_lambda_decay;
      Eigen::MatrixXd diffused_dfdV = sobolev_diffuse_gradients(dfdV, *tmp_solver.hullMesh, sobolev_lambda, sobolev_p);
      
      double raw_norm = tinyAD_flatten(dfdV).norm(); // use eigen flaten
      double diffused_norm = tinyAD_flatten(diffused_dfdV).norm();
      
      // std::cout << "energies: " << energies[0][iter] << " " << energies[1][iter] << " " << energies[2][iter] << "\n";
      curr_hull_psmesh->addVertexVectorQuantity(" sobolev diffused grads", diffused_dfdV*(-1.))->setEnabled(true);
      dfdV = diffused_dfdV;
    }
    //DEBUG
    polyscope::frameTick();
    // polyscope::screenshot(false);
    // polyscope::show();

    printf("line search\n");
    double opt_step_size = hull_update_line_search(dfdV, hull_positions, G_vec, fair_sides_count, 
                                                   dice_energy_step, dice_search_decay, frozen_G, 400);
    std::cout << ANSI_FG_RED << "  line search done!\n" << opt_step_size << ANSI_RESET << "\n";
    hull_positions = hull_positions - opt_step_size * dfdV;
    tmp_solver = Forward3DSolver(hull_positions, G_vec);
    // IMPORTANT step; tmo solver's conv hull will shuffle the order of vertices
    hull_positions = vertex_data_to_matrix(tmp_solver.hullGeometry->inputVertexPositions);
    tmp_solver.initialize_pre_computes();
  }
  // use boundary builder for visuals
  BoundaryBuilder tmp_bnd_builder(&tmp_solver);
  tmp_bnd_builder.build_boundary_normals();
  update_visuals(&tmp_solver, &tmp_bnd_builder);
  tmp_bnd_builder.print_area_of_boundary_loops();
}


// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
  if (ImGui::BeginCombo("##combo1", input_name.c_str())){ // The second parameter is the label previewed before opening the combo.
      for (std::string tmp_str: all_input_names){ // This enables not having to have a const char* arr[]. Or maybe I'm just a noob.
          bool is_selected = (input_name == tmp_str.c_str()); // You can store your selection however you want, outside or inside your objects
          if (ImGui::Selectable(tmp_str.c_str(), is_selected)){ // selected smth
              input_name = tmp_str;
              initialize_state(input_name);
          }
          if (is_selected)
              ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
      }
      ImGui::EndCombo();
  }
  
}


int main(int argc, char **argv) {
  // build mesh
  vis_utils = VisualUtils();
  polyscope::init();

  initialize_state(input_name);
  // Initialize polyscope
  // Set the callback function
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  polyscope::view::upDir = polyscope::view::UpDir::NegZUp;
  // init_convex_shape_to_fill(all_polygons_current_item2);
  // deformationSolver = new DeformationSolver(mesh, geometry, convex_to_fill_mesh, convex_to_fill_geometry);
  polyscope::state::userCallback = myCallback;

  polyscope::show();

  return EXIT_SUCCESS;
}
