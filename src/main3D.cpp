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
#include "visual_utils.h"
#include "markov_model.h"
#include "boundary_tools.h"

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
polyscope::SurfaceMesh *psInputMesh, *psHullMesh, *dummy_psMesh1, 
                       *dummy_psMesh2, *dummy_psMesh3, *dummy_psMesh4,
                       *dummy_forward_vis,
                       *coloredPsMesh, *gm_sphere_mesh, *height_function_mesh;

polyscope::PointCloud *curr_state_pt, *curr_g_vec_gm_pt,
                      *raster_pc;

polyscope::SurfaceGraphQuantity* curr_state_segment;

float pt_cloud_radi_scale = 0.01,
      curve_radi_scale = 0.002,
      G_r = 0.158,
      G_theta = 0.,
      G_phi = 0.,
      g_vec_theta = 0.,
      g_vec_phi = 0.;


bool gm_is_drawn = false,
     height_function_extras = false;
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

// snail trail stuff
bool show_snail_trail = true;
polyscope::PointCloud *snail_trail_pc;
polyscope::SurfaceMesh *dummy_ps_mesh_for_snail_trail;
Vector3 old_g_vec, new_g_vec;
int snail_trail_dummy_counter = 0;

// Markov Chain stuff
RollingMarkovModel *markov_model;
polyscope::PointCloud *sudo_faces_pc;

// boundary stuff
BoundaryBuilder *boundary_builder;
polyscope::PointCloud *boundary_normals_pc;
polyscope::SurfaceMesh *dummy_psMesh_for_regions, *dummy_psMesh_for_height_surface;
bool draw_boundary_patches = false;
bool test_guess = true;
polyscope::PointCloud *test_pc;

// example choice
std::vector<std::string> all_polyhedra_items = {std::string("cube"), std::string("tet"), std::string("sliced tet"), std::string("Conway spiral 4"), std::string("oloid"), std::string("gomboc"), std::string("bunny"), std::string("bunnylp")};
std::string all_polygons_current_item = "tet";
static const char* all_polygons_current_item_c_str = "tet";



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
  // just draw the sphere next to the main surface
  // std::vector<Vector3> sphere_pos = {vis_utils.center};
  // gauss_map_pc = polyscope::registerPointCloud("Gauss Map", sphere_pos);
  // gauss_map_pc->setPointColor({0.74,0.7,0.9});
  // gauss_map_pc->setPointRadius(gm_radi, false);
  // gauss_map_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);
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
}


void visualize_vertex_probabilities(){
  forwardSolver->compute_vertex_probabilities();
  polyscope::PointCloud* vertex_prob_cloud = polyscope::registerPointCloud("vertex probs", 
            forwardSolver->hullGeometry->inputVertexPositions);
    // set some options
  vertex_prob_cloud->addScalarQuantity("initial prob", forwardSolver->vertex_probabilities.toVector());
  vertex_prob_cloud->setPointRenderMode(polyscope::PointRenderMode::Sphere);
}

void generate_polyhedron_example(std::string poly_str, bool triangulate = false){
    // readManifoldSurfaceMesh()
  std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(poly_str);
  mesh = mesh_ptr.release();
  geometry = geometry_ptr.release();
  preprocess_mesh(mesh, geometry, triangulate );//|| std::strcmp(poly_str.c_str(), "gomboc") == 0
}


// initialize Markov chain and do the pre-computes
void initiate_markov_model(){
  markov_model = new RollingMarkovModel(forwardSolver);
  markov_model->initialize_pre_computes();
  markov_model->split_chain_edges_and_build_probability_pairs();
  markov_model->print_prob_pairs();
}


// visualize 
void visualize_sudo_faces(){
  std::vector<Vector3> sf_positions;
  size_t sf_count = 0;
  for (Halfedge he: markov_model->mesh->halfedges()){
    SudoFace *curr_sf = markov_model->root_sudo_face[he];
    if (curr_sf != nullptr){
      while (curr_sf->next_sudo_face != curr_sf){
        sf_positions.push_back(curr_sf->normal + vis_utils.center);
        sf_count++;
        curr_sf = curr_sf->next_sudo_face;
      }
    }
  }
  polyscope::PointCloud* sudo_faces_pc = polyscope::registerPointCloud("sudo faces", sf_positions);
  sudo_faces_pc->setEnabled(true);
  sudo_faces_pc->setPointRadius(vis_utils.gm_pt_radi * 1.0, false);
  sudo_faces_pc->setPointColor({0.,1.,1.});
}


// color input polyhedra with default
void color_faces_with_default(){
  FaceData<Vector3> face_colors(*forwardSolver->hullMesh, default_face_color);
  polyscope::SurfaceFaceColorQuantity *fColor = psHullMesh->addFaceColorQuantity("state vis color", face_colors);
  fColor->setEnabled(true);
}


// show current g vector
void visualize_g_vec(){
  std::vector<Vector3> the_g_vec = {forwardSolver->curr_g_vec};
  polyscope::PointCloudVectorQuantity *psG_vec = vis_utils.psG->addVectorQuantity("g_vec", the_g_vec);
  psG_vec->setEnabled(true);
  psG_vec->setVectorRadius(curve_radi_scale * 5., false);
  psG_vec->setVectorLengthScale(1., false);
  psG_vec->setVectorColor({0.8,0.1,0.1});
}


// visualize current touching element: Vertex/Edge/Face
void visualize_contact(){
  if (polyscope::hasPointCloud("current Vertex")) polyscope::removePointCloud("current Vertex");
  dummy_forward_vis->removeQuantity("current contact edge"); // has a built-in existance checker
  if (forwardSolver->curr_f.getIndex() == INVALID_IND) color_faces_with_default();
  // add the other two
  if (forwardSolver->curr_v.getIndex() != INVALID_IND) {
    printf("at Vertex vis\n");
    // show contact vertex on polyhedra
    std::vector<Vector3> curr_state_pos = {forwardSolver->hullGeometry->inputVertexPositions[forwardSolver->curr_v]}; // first and second should be the same since we just initialized.
    curr_state_pt = polyscope::registerPointCloud("current Vertex", curr_state_pos);
    curr_state_pt->setEnabled(true);
    curr_state_pt->setPointRadius(pt_cloud_radi_scale * 5., false);
    
    // show gravity vec on Gauss Map
    std::vector<Vector3> curr_g_vec_pos = {forwardSolver->curr_g_vec + vis_utils.center}; // first and second should be the same since we just initialized.
    curr_g_vec_gm_pt = polyscope::registerPointCloud("current g vec", curr_g_vec_pos);
    curr_g_vec_gm_pt->setEnabled(true);
    curr_g_vec_gm_pt->setPointRadius(vis_utils.gm_pt_radi * 5., false);
    curr_g_vec_gm_pt->setPointColor({0.,0.,0.});
  }
  else if (forwardSolver->curr_e.getIndex() != INVALID_IND){
    printf("at Edge vis\n");
    Vertex v1 = forwardSolver->curr_e.firstVertex(),
           v2 = forwardSolver->curr_e.secondVertex();
    Vector3 p1 = forwardSolver->hullGeometry->inputVertexPositions[v1],
            p2 = forwardSolver->hullGeometry->inputVertexPositions[v2];
    std::vector<std::array<size_t, 2>> edgeInds;
    std::vector<Vector3> positions;
    edgeInds.push_back({0, 1});
    positions.push_back(p1); positions.push_back(p2);
    curr_state_segment =  dummy_forward_vis->addSurfaceGraphQuantity("current contact edge", positions, edgeInds);
    curr_state_segment->setRadius(curve_radi_scale*3., false);
    curr_state_segment->setColor({0., 0., 1.});
    curr_state_segment->setEnabled(true);
  }
  else if (forwardSolver->curr_f.getIndex() != INVALID_IND){
    face_colors = FaceData<Vector3>(*forwardSolver->hullMesh, default_face_color);
    face_colors[forwardSolver->curr_f] = curr_face_color;
    polyscope::SurfaceFaceColorQuantity *fColor = psHullMesh->addFaceColorQuantity("state vis color", face_colors);
    fColor->setEnabled(true);
  }
  else {
    polyscope::warning("all elements are Invalid??\n");
  }
  // TODO: maybe add Snail trail?
}


void initialize_state_vis(){
  // for later single segment curve addition (for stable edges)
  std::vector<std::vector<size_t>> dummy_face{{0,0,0}};
  std::vector<Vector3> dummy_pos = {Vector3({0.,0.,0.})};
  dummy_forward_vis = polyscope::registerSurfaceMesh("state check mesh", dummy_pos, dummy_face); // nothing matters in this line
  // will add curves to this later
  vis_utils.draw_G();
  visualize_contact();
  visualize_g_vec();
  color_faces_with_default();
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
  auto t1 = clock();
  vis_utils.draw_G();
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
  
  if (ImGui::SliderFloat("G radi scale", &G_r, -1., 1.)||
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
    G = find_center_of_mass(*forwardSolver->hullMesh, *forwardSolver->hullGeometry).first;
    update_solver_and_boundaries();
    update_visuals_with_G();
  }
  if (ImGui::Button("Show Gauss Map")){///face_normal_vertex_gm_radi
    visualize_gauss_map();
    // Debugging stuff
    // geometry->requireFaceNormals();
    // FaceData<Vector3> selected_normals(*mesh, Vector3::zero());
    // selected_normals[mesh->face(566)] = geometry->faceNormal(mesh->face(566));
    // selected_normals[mesh->face(3975)] = geometry->faceNormal(mesh->face(3975));
    // psInputMesh->addFaceVectorQuantity("selected normals ", selected_normals)->setEnabled(true);
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
  if (ImGui::Checkbox("height function extras", &height_function_extras));
  
  // my simulation 
  if (ImGui::Button("initialize g_vec") ||
    ImGui::SliderFloat("initial g_vec theta", &g_vec_theta, 0., 2*PI)||
    ImGui::SliderFloat("initial g_vec phi", &g_vec_phi, 0., 2*PI)) {
    initial_g_vec = {cos(g_vec_phi)*sin(g_vec_theta), cos(g_vec_phi)*cos(g_vec_theta), sin(g_vec_phi)};
    // initial_g_vec = geometry->faceNormal(mesh->face(1942));
    forwardSolver->find_contact(initial_g_vec);
    initialize_state_vis();
    // snail trail stuff
    snail_trail_dummy_counter = 0;
    std::vector<std::vector<size_t>> dummy_face{{0,0,0}};
    dummy_ps_mesh_for_snail_trail = polyscope::registerSurfaceMesh("dummy mesh snail trail", geometry->inputVertexPositions, dummy_face);
    old_g_vec = forwardSolver->curr_g_vec;
  }
  if (ImGui::Button("next state")){
    forwardSolver->next_state(true);
    new_g_vec = forwardSolver->curr_g_vec;
    visualize_g_vec();
    visualize_contact();
    if(show_snail_trail && norm(old_g_vec-new_g_vec) != 0.){ // proly dont have to use tol
      draw_arc_on_sphere(old_g_vec, new_g_vec, vis_utils.center,
                         vis_utils.gm_radi, vis_utils.arcs_seg_count, 200 + snail_trail_dummy_counter, 
                         dummy_ps_mesh_for_snail_trail, 3., {0.9,0.8,0.1});
      snail_trail_dummy_counter++;
    }
    old_g_vec = new_g_vec;
  }
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
  if (ImGui::Button("show sudo faces")){
    initiate_markov_model();
    visualize_sudo_faces();
  }
  if (ImGui::Button("matrix debugg")){
    markov_model->build_transition_matrix();
    markov_model->check_transition_matrix();
  }
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
