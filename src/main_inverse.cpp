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
#include "polyscope/point_cloud.h"
#include "args/args.hxx"
#include "imgui.h"

#include "coloring.h"
#include "forward3D.h"
#include "mesh_factory.h"
// #include "geometry_utils.h" //; in fwd solver 
#include "visual_utils.h"
#include "markov_model.h"
#include "inv_design.h"
// #include "igl/arap.h"

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
bool draw_boundary_patches = true;
bool test_guess = true;
bool compute_global_G_effect = true,
     structured_opt = true,
     always_update_structure = true,
     with_hull_projection = false,
     use_initial_Ls = false,
     use_initial_geometry = false,
     first_time = false;

polyscope::PointCloud *test_pc;


// optimization stuff
InverseSolver* inverseSolver;
float step_size = 0.01,
      step_size2 = 0.01,
      step_size3 = 0.01;
float stable_normal_update_thresh = -1;
int ARAP_max_iters = 10;

// example choice
std::vector<std::string> all_polyhedra_items = {std::string("tet"), std::string("tet2"), std::string("cube"), std::string("tilted cube"), std::string("sliced tet"), std::string("Conway spiral 4"), std::string("D100"), std::string("D120"),std::string("oloid"), std::string("small_bunny"), std::string("bunnylp"), std::string("kitten"), std::string("double-torus"), std::string("soccerball"), std::string("cowhead"), std::string("bunny"), std::string("gomboc"), std::string("equi-180")};
std::string all_polygons_current_item = "bunnylp";
static const char* all_polygons_current_item_c_str = "bunnylp";


void draw_stable_patches_on_gauss_map(bool on_height_surface = false){
  if (!forwardSolver->updated)
    forwardSolver->initialize_pre_computes();
  // std::vector<Vector3> boundary_normals;
  if (draw_boundary_patches){
    auto net_pair = build_and_draw_stable_patches_on_gauss_map(boundary_builder, dummy_psMesh_for_regions,
                                               vis_utils.center, vis_utils.gm_radi, vis_utils.arcs_seg_count, 
                                               on_height_surface);
    if (test_guess)
      vis_utils.draw_guess_pc(net_pair.first, net_pair.second);
    else {
      if (polyscope::hasPointCloud("test point cloud")) polyscope::removePointCloud("test point cloud");
      if (polyscope::hasSurfaceMesh("dummy mesh for stable regions on polyhedra")) polyscope::removeSurfaceMesh("dummy mesh for stable regions on polyhedra");
    }
  }
  else
    if (polyscope::hasSurfaceMesh("dummy mesh for GM patch arcs")) polyscope::removeSurfaceMesh("dummy mesh for GM patch arcs");
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
  psHullMesh = polyscope::registerSurfaceMesh(
    "hull mesh",
    forwardSolver->hullGeometry->inputVertexPositions, forwardSolver->hullMesh->getFaceVertexList(),
    polyscopePermutations(*forwardSolver->hullMesh));
  vis_utils.draw_G();
  // psInputMesh->addFaceVectorQuantity("normals", geometry->faceNormals);
}


void generate_polyhedron_example(std::string poly_str, bool triangulate = false){
    // readManifoldSurfaceMesh()
  std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(poly_str);
  mesh = mesh_ptr.release();
  geometry = geometry_ptr.release();
  preprocess_mesh(mesh, geometry, triangulate || std::strcmp(poly_str.c_str(), "gomboc") == 0);
  first_time = true;
}


// color input polyhedra with default
void color_faces_with_default(){
  FaceData<Vector3> face_colors(*forwardSolver->hullMesh, default_face_color);
  polyscope::SurfaceFaceColorQuantity *fColor = psInputMesh->addFaceColorQuantity("state vis color", face_colors);
  fColor->setEnabled(true);
}


void color_faces(){
  if (recolor_faces){
    // printf("hull faces: %d\n", forwardSolver->hullMesh->nFaces());
    face_colors = generate_random_colors(forwardSolver->hullMesh);
    for (Face f: forwardSolver->hullMesh->faces()){
      if (!forwardSolver->face_is_stable(f))
        face_colors[f] = default_face_color;
    }
    recolor_faces = false;
  }
}

void update_solver_and_boundaries(){
  auto t1 = clock();
  forwardSolver->set_G(G);
  // printf("forward precomputes \n");
  forwardSolver->initialize_pre_computes();
  printf("building boundary normals \n");
  boundary_builder->build_boundary_normals();
  boundary_builder->print_area_of_boundary_loops();
  printf("precomputes took %d\n", clock() - t1);
}


void update_visuals_with_G(){
  
  vis_utils.draw_G();
  if (gm_is_drawn)
    vis_utils.plot_height_function();
  // stuff on the polyhedra
  auto t1 = clock();
  // visualize_edge_stability();
  // visualize_face_stability();
  // visualize_stable_vertices();
  
  // coloring stuff
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
  printf("GM visuals took %d\n", clock() - t1);
  t1 = clock();
  // initialize_boundary_builder();
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
  inverseSolver->find_d_pf_d_Gs();
  Vector3 G_grad = inverseSolver->find_total_g_grad();
  G = G + step_size * G_grad;
  update_solver_and_boundaries();
  update_visuals_with_G();
}


void take_opt_vertices_step(){
  inverseSolver->find_d_pf_dvs();
  VertexData<Vector3> total_vertex_grads = inverseSolver->find_total_vertex_grads();
  VertexData<Vector3> new_poses(*forwardSolver->hullMesh);
  for (Vertex v: forwardSolver->hullMesh->vertices()){
    forwardSolver->hullGeometry->inputVertexPositions[v] += step_size2 * total_vertex_grads[v];
  }
  polyscope::getSurfaceMesh("hull mesh")->updateVertexPositions(forwardSolver->hullGeometry->inputVertexPositions);
  update_solver_and_boundaries();
  update_visuals_with_G();
}


void take_uni_mass_opt_vertices_step(){

  printf("finding uniform mass\n");
  forwardSolver->set_uniform_G();
  printf("finding uniform mass\n");
  forwardSolver->initialize_pre_computes();
  printf("finding derivatives\n");
  inverseSolver->find_uni_mass_d_pf_dv();
  printf(" uni mass derivatives found!\n");
  if (structured_opt){ // if updating the flow structure
    boundary_builder->build_boundary_normals();
    boundary_builder->print_area_of_boundary_loops();
    if (first_time || always_update_structure){
      inverseSolver->flow_structure = forwardSolver->face_last_face;
      if (first_time) stable_normal_update_thresh = -1.;
      else stable_normal_update_thresh = 0.1;
      first_time = false;
    }
  }
  printf(" finding total per vertex derivatives\n");
  VertexData<Vector3> total_uni_mass_vertex_grads = inverseSolver->find_uni_mass_total_vertex_grads(structured_opt, stable_normal_update_thresh);
  printf(" total per vertex derivatives found!\n");

  // subdivide for aggressive updates
  // inverseSolver->subdivide_for_aggressive_updates(total_uni_mass_vertex_grads);

  // Update original mesh vertices; hull will be updated internally
  // polyscope::registerSurfaceMesh("pre-step hull", forwardSolver->hullGeometry->inputVertexPositions,
  //                                              forwardSolver->hullMesh->getFaceVertexList());
  polyscope::registerSurfaceMesh("pre-step mesh", forwardSolver->inputGeometry->inputVertexPositions,
                                                  forwardSolver->inputMesh->getFaceVertexList())->setEnabled(false);
  
  // PC
  polyscope::PointCloud* pre_pc = polyscope::registerPointCloud("pre-step hull points", forwardSolver->hullGeometry->inputVertexPositions);
  pre_pc->setPointRadius(vis_utils.gm_pt_radi * 1., false);
  pre_pc->setPointColor(glm::vec3({0.8,0.1, 0.1}));
  pre_pc->setEnabled(false);
  
  printf("update interior indexing\n");
  forwardSolver->update_hull_points_correspondence(forwardSolver->hullGeometry->inputVertexPositions 
                                                   + total_uni_mass_vertex_grads * step_size3, 
                                                   inverseSolver->initial_geometry->inputVertexPositions);
  
  VertexData<Vector3> trivial_updates;
  trivial_updates = inverseSolver->trivial_update_positions(total_uni_mass_vertex_grads * step_size3);
  auto new_positions = forwardSolver->inputGeometry->inputVertexPositions + trivial_updates;
  polyscope::registerSurfaceMesh("pre-ARAP mesh", new_positions,
                                               forwardSolver->inputMesh->getFaceVertexList())->setEnabled(false);
  // update indexing for interior vertices


  printf("updating interior positions\n");
  VertexData<Vector3> updates;
  if (forwardSolver->hullMesh->nVertices() == forwardSolver->inputMesh->nVertices()) // no interior vertices
    updates = trivial_updates;
  else {
    // updates = inverseSolver->laplace_update_positions(new_positions);
    updates = inverseSolver->ARAP_update_positions(new_positions);
    // updates = inverseSolver->greedy_update_positions(total_uni_mass_vertex_grads * step_size3);
    // updates = inverseSolver->diffusive_update_positions(total_uni_mass_vertex_grads * step_size3);
  }
  forwardSolver->inputGeometry->inputVertexPositions += updates;
  printf("updates for interior found!\n");
  // update hull with new positions
  forwardSolver->update_convex_hull(with_hull_projection);
  
  polyscope::registerSurfaceMesh("post-ARAP mesh", forwardSolver->inputGeometry->inputVertexPositions,
                                               forwardSolver->inputMesh->getFaceVertexList());
  polyscope::registerSurfaceMesh("hull mesh", forwardSolver->hullGeometry->inputVertexPositions,
                                              forwardSolver->hullMesh->getFaceVertexList())->setEnabled(false);
  // PC
  polyscope::PointCloud* post_pc = polyscope::registerPointCloud("post-step hull points", forwardSolver->hullGeometry->inputVertexPositions);
  pre_pc->setPointRadius(vis_utils.gm_pt_radi * 1., false);
  pre_pc->setPointColor(glm::vec3({0.1,0.1, 0.8}));
  forwardSolver->set_uniform_G();
  G = forwardSolver->get_G();
  update_solver_and_boundaries();
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

              //
              visualize_gauss_map();//
              G = find_center_of_mass(*forwardSolver->inputMesh, *forwardSolver->inputGeometry).first;
              update_solver_and_boundaries();
              update_visuals_with_G();
          }
          if (is_selected)
              ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
      }
      ImGui::EndCombo();
  }

  if (ImGui::SliderFloat("vertex radi scale", &pt_cloud_radi_scale, 0., 1.)){
    vis_utils.draw_G();
  }
  
  if (ImGui::SliderFloat("G radi scale", &G_r, -1., 1.)||
      ImGui::SliderFloat("G theta", &G_theta, 0., 2*PI)||
      ImGui::SliderFloat("G phi", &G_phi, 0., 2*PI)) {
    G = {cos(G_phi)*sin(G_theta)*G_r, cos(G_phi)*cos(G_theta)*G_r, sin(G_phi)*G_r};
    if (G_is_inside(*forwardSolver->hullMesh, *forwardSolver->hullGeometry, G)){
      update_solver_and_boundaries();
      update_visuals_with_G();
      if (real_time_raster){
        color_faces();
        build_raster_image();
      }
    }
  }
  if (ImGui::Button("uniform mass G")){
    G = find_center_of_mass(*forwardSolver->inputMesh, *forwardSolver->inputGeometry).first;
    update_solver_and_boundaries();
    update_visuals_with_G();
  }
  if (ImGui::Checkbox("colored arcs (slows)", &color_arcs)) vis_utils.color_arcs = color_arcs;
  if (ImGui::Checkbox("show hidden stable vertex normals", &show_hidden_stable_vertex_normals)) {
    vis_utils.show_hidden_stable_vertex_normals = show_hidden_stable_vertex_normals;
    vis_utils.draw_stable_vertices_on_gauss_map();
  }
  // my simulation 
  if (ImGui::SliderInt("sample count", &sample_count, 1000, 1000000));
  if (ImGui::Button("build raster image")){
    color_faces();
    vis_utils.visualize_colored_polyhedra(face_colors);
    build_raster_image();
  }
  // if(ImGui::Button("recolor faces")){
  //   recolor_faces = true;// so bad lol
  //   color_faces();
  //   vis_utils.visualize_colored_polyhedra(face_colors);
  // }
  if (ImGui::Checkbox("draw artificial R3 boundaries", &test_guess)) draw_stable_patches_on_gauss_map();
  if (ImGui::Button("FD (move G)")) {
    inverseSolver->find_d_pf_d_Gs(true);
  }
  if (ImGui::Button("FD (move vertices)")) {
    inverseSolver->find_d_pf_dvs(true);
  }
  if (ImGui::Button("joint gradient (FD)")) {
    forwardSolver->set_uniform_G();
    forwardSolver->initialize_pre_computes();
    inverseSolver->find_uni_mass_d_pf_dv(true);
  }
  if (ImGui::Button("take fair step (joint gradient)")) {
    take_uni_mass_opt_vertices_step();
  }
  if (ImGui::SliderFloat("step size 3", &step_size3, 0., 1.0));
  if (ImGui::Checkbox("compute global G effect", &compute_global_G_effect)) inverseSolver->compute_global_G_effect = compute_global_G_effect;
  if (ImGui::Checkbox("structured opt", &structured_opt));
  if (ImGui::Checkbox("update structured at every step", &always_update_structure));
  if (ImGui::Checkbox("decrease convex hull points", &with_hull_projection));
  if (ImGui::Checkbox("use initial geomtery", &use_initial_geometry)) {
    inverseSolver->use_old_geometry = use_initial_geometry;
  }
  if (ImGui::SliderInt("ARAP max iters", &ARAP_max_iters, 1, 30)) inverseSolver->arap_max_iter = ARAP_max_iters;
  if (ImGui::SliderFloat("stable normal update thresh", &stable_normal_update_thresh, 0., 4.0));

}


int main(int argc, char **argv) {
  // build mesh
  vis_utils = VisualUtils();
  generate_polyhedron_example(all_polygons_current_item);
  G = {0.,0,0};
  update_solver();
  // convex_hull(forwardSolver->hullGeometry->inputVertexPositions);
  // build the solver

  // Initialize polyscope
  polyscope::init();
  color_faces();
  vis_utils.visualize_colored_polyhedra(face_colors);
  visualize_gauss_map();
  // Set the callback function
  polyscope::state::userCallback = myCallback;
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  polyscope::view::upDir = polyscope::view::UpDir::NegZUp;
  
  // Give control to the polyscope gui
  polyscope::show();
  

  return EXIT_SUCCESS;
}
