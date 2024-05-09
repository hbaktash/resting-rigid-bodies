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
std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr, cv_mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr, cv_geometry_ptr;
ManifoldSurfaceMesh* mesh, *convex_to_fill_mesh;
VertexPositionGeometry* geometry, *convex_to_fill_geometry;
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
bool draw_boundary_patches = false;
bool test_guess = true;
bool compute_global_G_effect = true,
     deform_after = true,
     frozen_G = false,
     structured_opt = true,
     one_time_CP_assignment = false,
     always_update_structure = true,
     with_hull_projection = false,
     first_time = false;
polyscope::PointCloud *test_pc;

int fair_sides_count = 6; // for optimization


// optimization stuff
InverseSolver* inverseSolver;
float step_size = 0.01,
      step_size2 = 0.01,
      step_size3 = 0.167;
float stable_normal_update_thresh = -1;
int ARAP_max_iters = 10;

// deformation
bool animate = false,
     v2_dice_animate = false;
float scale_for_feasi = 2.;
float membrane_lambda = 0.2,
      CP_lambda_exp = 1.,
      CP_mu = 1.1,
      refinement_CP_threshold = 0.3,
      barrier_lambda_exp = 1.,
      barrier_mu = 0.8;
int filling_max_iter = 200;
int hull_opt_steps = 50;
// example choice
std::vector<std::string> all_polyhedra_items = {std::string("tet"), std::string("tet2"), std::string("cube"), std::string("tilted cube"), std::string("sliced tet"), std::string("Conway spiral 4"), std::string("oloid"), std::string("fox"), std::string("small_bunny"), std::string("bunnylp"), std::string("kitten"), std::string("double-torus"), std::string("soccerball"), std::string("cowhead"), std::string("bunny"), std::string("gomboc"), std::string("mark_gomboc")};
std::string all_polygons_current_item = "sliced tet",
            all_polygons_current_item2 = "tet";
static const char* all_polygons_current_item_c_str = "bunnylp";

// oris
std::vector<std::string> all_stable_face_items = {"nothing yet"};
std::string all_stable_face_current_item = "nothing yet";
static const char* all_stable_face_current_item_c_str = "nothing yet";


void draw_stable_patches_on_gauss_map(bool on_height_surface = false, 
                                      BoundaryBuilder *bnd_builder = boundary_builder,
                                      Forward3DSolver *tmp_solver = forwardSolver){
  if (!forwardSolver->updated)
    tmp_solver->initialize_pre_computes();
  // std::vector<Vector3> boundary_normals;
  if (draw_boundary_patches){
    auto net_pair = build_and_draw_stable_patches_on_gauss_map(bnd_builder, dummy_psMesh_for_regions,
                                               vis_utils.center, vis_utils.gm_radi, vis_utils.arcs_seg_count, 
                                               on_height_surface);
    if (test_guess)
      vis_utils.draw_guess_pc(net_pair.first, net_pair.second);
    else {
      if (polyscope::hasPointCloud("test point cloud")) polyscope::removePointCloud("test point cloud");
      if (polyscope::hasCurveNetwork("stable regions on polyhedra")) polyscope::removeCurveNetwork("stable regions on polyhedra");
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
  vis_utils.draw_G();
  // psInputMesh->addFaceVectorQuantity("normals", geometry->faceNormals);
}

void update_solver(){
  //assuming convex input here
  forwardSolver = new Forward3DSolver(mesh, geometry, G);
  forwardSolver->initialize_pre_computes();
  initialize_boundary_builder();
  inverseSolver = new InverseSolver(boundary_builder);
  vis_utils.forwardSolver = forwardSolver;
}

void generate_polyhedron_example(std::string poly_str, bool triangulate = false){
    // readManifoldSurfaceMesh()
  std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(poly_str);
  mesh = mesh_ptr.release();
  geometry = geometry_ptr.release();
  preprocess_mesh(mesh, geometry, triangulate || std::strcmp(poly_str.c_str(), "gomboc") == 0);
  first_time = true;
}


void color_faces(Forward3DSolver *fwd_solver){
  if (recolor_faces){
    // printf("hull faces: %d\n", forwardSolver->hullMesh->nFaces());
    face_colors = generate_random_colors(fwd_solver->hullMesh);
    for (Face f: fwd_solver->hullMesh->faces()){
      if (!fwd_solver->face_is_stable(f))
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
  // boundary_builder->print_area_of_boundary_loops();
  // printf("precomputes took %d\n", clock() - t1);
}


void update_visuals_with_G(Forward3DSolver *tmp_solver = nullptr, BoundaryBuilder *bnd_builder = boundary_builder){
  if (tmp_solver != nullptr)
    vis_utils.forwardSolver = tmp_solver;
  
  vis_utils.draw_G();
  if (gm_is_drawn)
    vis_utils.plot_height_function();
  // stuff on the polyhedra
  auto t1 = clock();
  
  // coloring stuff
  recolor_faces = true;// so bad lol
  color_faces(tmp_solver);
  vis_utils.visualize_colored_polyhedra(face_colors);
  // // Gauss map stuff
  if(color_arcs){
    vis_utils.draw_edge_arcs_on_gauss_map();
  }
  vis_utils.draw_stable_vertices_on_gauss_map();
  vis_utils.show_edge_equilibria_on_gauss_map();
  vis_utils.draw_stable_face_normals_on_gauss_map();
  // printf("GM visuals took %d\n", clock() - t1);
  if (draw_boundary_patches)
    draw_stable_patches_on_gauss_map(gm_is_drawn, bnd_builder);
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
              polyscope::removeAllStructures();
              all_polygons_current_item = tmp_str;
              generate_polyhedron_example(all_polygons_current_item);
              update_solver();
              init_visuals();
              recolor_faces = true;
              color_faces(forwardSolver);

              // //
              visualize_gauss_map();//
              G = find_center_of_mass(*forwardSolver->inputMesh, *forwardSolver->inputGeometry).first;
              update_solver_and_boundaries();
              boundary_builder->print_area_of_boundary_loops();
              printf("here1\n");
              update_visuals_with_G(forwardSolver, boundary_builder);
              printf("here2\n");
              if (polyscope::hasSurfaceMesh("fillable hull")) polyscope::getSurfaceMesh("fillable hull")->setEnabled(false);
              if (polyscope::hasSurfaceMesh("temp sol")) polyscope::getSurfaceMesh("temp sol")->setEnabled(false);
          }
          if (is_selected)
              ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
      }
      ImGui::EndCombo();
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
    color_faces(forwardSolver);
    vis_utils.visualize_colored_polyhedra(face_colors);
    build_raster_image();
  }
  if (ImGui::Checkbox("draw artificial R3 boundaries", &test_guess)) draw_stable_patches_on_gauss_map();
  if (ImGui::Button("show all stable oris (sorted by prob)")){
    // vis_utils.visualize_all_stable_orientations();
    //
    forwardSolver->set_uniform_G();
    forwardSolver->initialize_pre_computes();
    BoundaryBuilder *bnd_builder = new BoundaryBuilder(forwardSolver);
    bnd_builder->build_boundary_normals();

    std::vector<std::pair<Face, double>> probs;
    for (Face f: forwardSolver->hullMesh->faces())
        if (bnd_builder->face_region_area[f] > 0)
            probs.push_back({f, bnd_builder->face_region_area[f]/(4.*PI)});
    std::sort(probs.begin(), probs.end(), [] (auto a, auto b) { return a.second > b.second; });
    all_stable_face_items = {};
    for (auto p: probs){
      Face f = p.first;
      double prob = p.second;
      all_stable_face_items.push_back("f " + std::to_string(f.getIndex()) + " prob: " + std::to_string(prob));
    }
  }
  if (ImGui::BeginCombo("##combo2", all_stable_face_current_item.c_str())){
      for (std::string tmp_str: all_stable_face_items){ 
          bool is_selected = (all_stable_face_current_item == tmp_str.c_str()); // You can store your selection however you want, outside or inside your objects
          if (ImGui::Selectable(tmp_str.c_str(), is_selected)){ // selected smth
              all_stable_face_current_item = tmp_str;
              
              size_t pos = 0;
              std::string token;
              pos = tmp_str.find(" ");
              token = tmp_str.substr(0, pos);
              tmp_str.erase(0, pos + 1);
              pos = tmp_str.find(" ");
              token = tmp_str.substr(0, pos);
              size_t f_ind = std::stoi(token);
              
              double floor_z = -1.;
              Vector3 floor_vec = Vector3({0,0,floor_z});
              VertexData<Vector3> rotated_poses(*forwardSolver->inputMesh);
              size_t i = 0;
              Face f = forwardSolver->hullMesh->face(f_ind);
              //   double prob = p.second;
              Vector3 f_normal = forwardSolver->hullGeometry->faceNormal(f);
              Vector3 rot_axis = cross(f_normal, floor_vec);
              double rot_angle = angle(f_normal, floor_vec);
              Vector3 offset = floor_vec - forwardSolver->hullGeometry->inputVertexPositions[f.halfedge().vertex()].rotateAround(rot_axis, rot_angle);
              for (Vertex v: forwardSolver->inputMesh->vertices()){
                rotated_poses[v] = forwardSolver->inputGeometry->inputVertexPositions[v].rotateAround(rot_axis, rot_angle) + offset;
              }
              double y_shift = 0.;
              Vector3 vis_shift({y_shift * 2., -2., 0.});
              auto tmp_ori_mesh = polyscope::registerSurfaceMesh("tmp orientation", -1. * rotated_poses + vis_shift, forwardSolver->inputMesh->getFaceVertexList());
              tmp_ori_mesh->setSurfaceColor({0.1, 0.4, 0.02});
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

  generate_polyhedron_example(all_polygons_current_item);
  G = {0.,0,0};
  update_solver();
  init_visuals();

  // Initialize polyscope
  // Set the callback function
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  polyscope::view::upDir = polyscope::view::UpDir::NegZUp;
  polyscope::state::userCallback = myCallback;

  polyscope::show();
  // Give control to the polyscope gui
  

  return EXIT_SUCCESS;
}
