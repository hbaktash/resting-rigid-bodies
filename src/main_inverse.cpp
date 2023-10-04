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

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
// #include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
// #include "bullet3/examples/BasicExample.h"
#include "args/args.hxx"
#include "imgui.h"

#include "coloring.h"
#include "forward3D.h"
#include "mesh_factory.h"
#include "geometry_utils.h"
#include "visual_utils.h"
#include "markov_model.h"
#include "inv_design.h"

// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Polyhedron_3.h>
// #include <CGAL/Surface_mesh.h>
// #include <CGAL/convex_hull_3.h>
// #include <vector>
// #include <fstream>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr;
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;
Vector3 G, // center of Mass
        initial_g_vec({0,-1,0}),
        default_face_color({0.99,0.99,0.99}),
        curr_face_color({0.1,0.87,0.1});

Forward3DSolver* forwardSolver;

VisualUtils vis_utils;

ManifoldSurfaceMesh* sphere_mesh;
VertexPositionGeometry* sphere_geometry;
        

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psInputMesh;

polyscope::PointCloud *raster_pc;

float pt_cloud_radi_scale = 0.1,
      G_r = 0.,
      G_theta = 0.,
      G_phi = 0.;


bool gm_is_drawn;
bool color_arcs = false,
     draw_unstable_edge_arcs = true,
     show_hidden_stable_vertex_normals = false;

// raster image stuff
int sample_count = 8000; 
FaceData<Vector3> face_colors;
bool recolor_faces = true,
     real_time_raster = false,
     draw_point_cloud = false,
     normalize_vecF = true;
glm::vec3 vecF_color({0.1, 0.1, 0.1});

// boundary stuff
BoundaryBuilder *boundary_builder;
polyscope::PointCloud *boundary_normals_pc;
polyscope::SurfaceMesh *dummy_psMesh_for_regions, *dummy_psMesh_for_height_surface;
bool draw_boundary_patches = false;
bool test_guess = true;
polyscope::PointCloud *test_pc;


// optimization stuff
InverseSolver* inverseSolver;
float step_size = 0.01;


// example choice
std::vector<std::string> all_polyhedra_items = {std::string("cube"), std::string("tilted cube"), std::string("tet"), std::string("sliced tet"), std::string("Conway spiral 4")};
std::string all_polygons_current_item = "tet";
static const char* all_polygons_current_item_c_str = "tet";



void draw_stable_patches_on_gauss_map(bool on_height_surface = false){
  if (!forwardSolver->updated)
    forwardSolver->initialize_pre_computes();
  std::vector<Vector3> boundary_normals;
  if (draw_boundary_patches){
    boundary_normals = build_and_draw_stable_patches_on_gauss_map(boundary_builder, dummy_psMesh_for_regions,
                                               vis_utils.center, vis_utils.gm_radi, vis_utils.arcs_seg_count, 
                                               on_height_surface);
  }
  if (test_guess)
    vis_utils.draw_guess_pc(boundary_normals);
}

//

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


// initialize boundary builder
void initialize_boundary_builder(){
  boundary_builder = new BoundaryBuilder(forwardSolver);
  boundary_builder->build_boundary_normals();
}

void update_solver(){
  //assuming convex input here
  forwardSolver = new Forward3DSolver(mesh, geometry, G);
  forwardSolver->initialize_pre_computes();
  initialize_boundary_builder();
  inverseSolver = new InverseSolver(boundary_builder);
  vis_utils.forwardSolver = forwardSolver;
  // Register the mesh with polyscope
  psInputMesh = polyscope::registerSurfaceMesh(
      "input mesh",
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));
  psInputMesh->setTransparency(0.75);
  vis_utils.draw_G();
}


// generate simple examples
void generate_polyhedron_example(std::string poly_str){
    // readManifoldSurfaceMesh()
    std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(poly_str);
    mesh = mesh_ptr.release();
    geometry = geometry_ptr.release();
}

// visualize boundary


// color input polyhedra with default
void color_faces_with_default(){
  FaceData<Vector3> face_colors(*forwardSolver->hullMesh, default_face_color);
  polyscope::SurfaceFaceColorQuantity *fColor = psInputMesh->addFaceColorQuantity("state vis color", face_colors);
  fColor->setEnabled(true);
}


void color_faces(){
  if (recolor_faces){
    face_colors = generate_random_colors(forwardSolver->hullMesh);
    for (Face f: forwardSolver->hullMesh->faces()){
      if (!forwardSolver->face_is_stable(f))
        face_colors[f] = default_face_color;
    }
    recolor_faces = false;
  }
}


void update_visuals_with_G(){
  forwardSolver->set_G(G);
  vis_utils.draw_G();
  auto t1 = clock();
  forwardSolver->initialize_pre_computes();
  printf("precomputes took %d\n", clock() - t1);
  if (gm_is_drawn)
    vis_utils.plot_height_function();
  // stuff on the polyhedra
  t1 = clock();
  // visualize_edge_stability();
  // visualize_face_stability();
  // visualize_stable_vertices();
  // // coloring stuff
  recolor_faces = true;// so bad lol
  color_faces();
  vis_utils.visualize_colored_polyhedra(face_colors);
  // // Gauss map stuff
  if(color_arcs){
    vis_utils.draw_edge_arcs_on_gauss_map();
  }
  vis_utils.draw_stable_vertices_on_gauss_map();
  vis_utils.show_edge_equilibria_on_gauss_map();
  vis_utils.draw_stable_face_normals_on_gauss_map();
  printf("visuals took %d\n", clock() - t1);
  t1 = clock();
  // initialize_boundary_builder();
  boundary_builder->build_boundary_normals();
  boundary_builder->print_area_of_boundary_loops();
  printf("boundary took %d\n", clock() - t1);
  t1 = clock();
  if (draw_boundary_patches)
    draw_stable_patches_on_gauss_map(gm_is_drawn);
  printf("boundary visuals took %d\n", clock() - t1);
  // initiate_markov_model();
  // visualize_sudo_faces();
}


// maybe move to another file, to merge with bullet sim; this and some other functions
// sample and raster; colors should be generated beforehand
void build_raster_image(){
  FaceData<std::vector<Vector3>> face_samples(*forwardSolver->hullMesh);
  FaceData<double> empirical_face_probabilities(*forwardSolver->hullMesh, 0.);
  int total_invalids = 0,
      valid_samples = 0;
  // need to re-iterate later with the same order
  FaceData<std::vector<Vector3>> initial_rolling_dirs(*forwardSolver->hullMesh);
  for (int i = 0; i < sample_count; i++){
    Vector3 random_g_vec = {randomReal(-1,1), randomReal(-1,1), randomReal(-1,1)};
    if (random_g_vec.norm() <= 1.){
      // if (i % 10000 == 0 && i > 0)
      //   printf("$$$ at sample %d\n", i);
      random_g_vec /= norm(random_g_vec);
      Face touching_face = forwardSolver->final_touching_face(random_g_vec);
      if (touching_face.getIndex() == INVALID_IND){
        total_invalids++;
        continue;
      }
      valid_samples++;
      empirical_face_probabilities[touching_face] += 1;
      if (draw_point_cloud){
        face_samples[touching_face].push_back(random_g_vec);
        if (normalize_vecF)
          initial_rolling_dirs[touching_face].push_back(forwardSolver->initial_roll_dir.normalize());
        else 
          initial_rolling_dirs[touching_face].push_back(forwardSolver->initial_roll_dir);
      }
    }
  }
  // printf(" ###### total invalid faces: %d  ######\n", total_invalids);
  
  std::vector<Vector3> raster_positions,
                       raster_colors;
  std::vector<Vector3> all_rolling_dirs;
  for (Face f: forwardSolver->hullMesh->faces()){
    empirical_face_probabilities[f] /= (double)valid_samples;
    if (draw_point_cloud) {
      std::vector<Vector3> tmp_points = face_samples[f];
      for (Vector3 tmp_p: tmp_points){
        raster_positions.push_back(tmp_p + vis_utils.center);
        raster_colors.push_back(face_colors[f]);
      }
      std::vector<Vector3> tmp_rolling_dirs = initial_rolling_dirs[f];
      for (Vector3 rolling_dir: tmp_rolling_dirs){
        all_rolling_dirs.push_back(rolling_dir);
      }
    }
  }
  if (draw_point_cloud) {
    raster_pc = polyscope::registerPointCloud("raster point cloud", raster_positions);
    raster_pc->setPointRadius(0.0078, false);
    polyscope::PointCloudColorQuantity* pc_col_quant = raster_pc->addColorQuantity("random color", raster_colors);
    pc_col_quant->setEnabled(true);
    polyscope::PointCloudVectorQuantity* raster_pc_vec_field = raster_pc->addVectorQuantity("init roll directions", all_rolling_dirs);
    raster_pc_vec_field->setEnabled(true);
    raster_pc_vec_field->setVectorLengthScale(0.05, false);
    raster_pc_vec_field->setVectorRadius(0.003, false);
    raster_pc_vec_field->setVectorColor(vecF_color);
  }
  std::cout << "empirical_face_probabilities: \n"<< empirical_face_probabilities.toVector().transpose() << "sum: " << empirical_face_probabilities.toVector().sum()<< "\n";
}


// inverse stuff
void take_opt_G_step(){
  inverseSolver->find_per_face_G_grads();
  Vector3 G_grad = inverseSolver->find_total_g_grad();
  G = G + step_size * G_grad;
  update_visuals_with_G();
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
  if (ImGui::BeginCombo("##combo1", all_polygons_current_item.c_str())){ // The second parameter is the label previewed before opening the combo.
      for (std::string tmp_str: all_polyhedra_items){ // This enables not having to have a const char* arr[]. Or maybe I'm just a noob.
          bool is_selected = (all_polygons_current_item == tmp_str.c_str()); // You can store your selection however you want, outside or inside your objects
          if (ImGui::Selectable(tmp_str.c_str(), is_selected)){ // selected smth
              all_polygons_current_item = tmp_str;
              generate_polyhedron_example(all_polygons_current_item);
              update_solver();
              recolor_faces = true;
              color_faces();
          }
          if (is_selected)
              ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
      }
      ImGui::EndCombo();
  }

  if (ImGui::SliderFloat("vertex radi scale", &pt_cloud_radi_scale, 0., 1.)){
    vis_utils.draw_G();
  }
  
  if (ImGui::SliderFloat("G radi scale", &G_r, -3., 3.)||
      ImGui::SliderFloat("G theta", &G_theta, 0., 2*PI)||
      ImGui::SliderFloat("G phi", &G_phi, 0., 2*PI)) {
    G = {cos(G_phi)*sin(G_theta)*G_r, cos(G_phi)*cos(G_theta)*G_r, sin(G_phi)*G_r};
    if (G_is_inside(*forwardSolver->hullMesh, *forwardSolver->hullGeometry, G)){
      update_visuals_with_G();
      if (real_time_raster){
        color_faces();
        build_raster_image();
      }
    }
  }
  if (ImGui::Button("uniform mass G")){
    G = find_center_of_mass(*mesh, *geometry);
    update_visuals_with_G();
  }
  if (ImGui::Button("Show Gauss Map")){///face_normal_vertex_gm_radi
    visualize_gauss_map();
  }
  if (ImGui::Checkbox("colored arcs (slows)", &color_arcs)) vis_utils.color_arcs = color_arcs;
  if (ImGui::Checkbox("draw unstable edge arcs", &draw_unstable_edge_arcs)) {
    vis_utils.draw_unstable_edge_arcs = draw_unstable_edge_arcs;
    vis_utils.draw_edge_arcs_on_gauss_map();
  }
  if (ImGui::Checkbox("show hidden stable vertex normals", &show_hidden_stable_vertex_normals)) {
    vis_utils.show_hidden_stable_vertex_normals = show_hidden_stable_vertex_normals;
    vis_utils.draw_stable_vertices_on_gauss_map();
  }
  if (ImGui::Button("Draw patches")){
    draw_stable_patches_on_gauss_map();
  }
  
  // my simulation 
  if (ImGui::SliderInt("sample count", &sample_count, 1000, 1000000));
  if (ImGui::Button("build raster image")){
    color_faces();
    vis_utils.visualize_colored_polyhedra(face_colors);
    build_raster_image();
  }
  if(ImGui::Button("recolor faces")){
    recolor_faces = true;// so bad lol
    color_faces();
    vis_utils.visualize_colored_polyhedra(face_colors);
  }
  if (ImGui::Checkbox("draw point cloud", &draw_point_cloud));
  if (ImGui::Checkbox("real time raster", &real_time_raster));
  if (ImGui::Checkbox("normalize vector field", &normalize_vecF));

  if (ImGui::Button("build and show region boundaries")) {
    initialize_boundary_builder();
    printf(" - pre patch draw!\n");
    draw_stable_patches_on_gauss_map();
  }
  if (ImGui::Checkbox("draw boundary patches", &draw_boundary_patches)) draw_stable_patches_on_gauss_map();
  if (ImGui::Checkbox("test guess", &test_guess)) draw_stable_patches_on_gauss_map();
  if (ImGui::Button("take fair step (move G)")) {
    take_opt_G_step();
  }
  if (ImGui::SliderFloat("step size", &step_size, 0., 0.10));
}


int main(int argc, char **argv) {
  // build mesh
  vis_utils = VisualUtils();
  generate_polyhedron_example(all_polygons_current_item);
  G = {0.,0,0};
  update_solver();
  // build the solver

  // Initialize polyscope
  polyscope::init();
  color_faces();
  vis_utils.visualize_colored_polyhedra(face_colors);
  // Set the callback function
  polyscope::state::userCallback = myCallback;
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  polyscope::view::upDir = polyscope::view::UpDir::NegZUp;
  
  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
