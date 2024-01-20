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
bool draw_boundary_patches = true;
bool test_guess = true;
bool compute_global_G_effect = true,
     line_search = false,
     deform_after = true,
     structured_opt = true,
     one_time_CP_assignment = true,
     always_update_structure = true,
     with_hull_projection = false,
     use_initial_geometry = false,
     first_time = false;

polyscope::PointCloud *test_pc;

int fair_sides_count = 6; // for optimization


// optimization stuff
InverseSolver* inverseSolver;
float step_size = 0.01,
      step_size2 = 0.01,
      step_size3 = 0.01;
float stable_normal_update_thresh = -1;
int ARAP_max_iters = 10;

// deformation
DeformationSolver *deformationSolver;
float CP_lambda_exp = 1.,
      CP_mu = 1.;
int filling_max_iter = 20;

// example choice
std::vector<std::string> all_polyhedra_items = {std::string("tet"), std::string("tet2"), std::string("cube"), std::string("tilted cube"), std::string("sliced tet"), std::string("Conway spiral 4"), std::string("oloid"), std::string("small_bunny"), std::string("bunnylp"), std::string("kitten"), std::string("double-torus"), std::string("soccerball"), std::string("cowhead"), std::string("bunny"), std::string("gomboc")};
std::string all_polygons_current_item = "sliced tet",
            all_polygons_current_item2 = "tet";
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

void init_visuals(){
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

void update_solver(){
  //assuming convex input here
  forwardSolver = new Forward3DSolver(mesh, geometry, G);
  forwardSolver->initialize_pre_computes();
  initialize_boundary_builder();
  inverseSolver = new InverseSolver(boundary_builder);
  vis_utils.forwardSolver = forwardSolver;
}


void init_convex_shape_to_fill(std::string poly_str, bool triangulate = false){
    // readManifoldSurfaceMesh()
  std::tie(cv_mesh_ptr, cv_geometry_ptr) = generate_polyhedra(poly_str);
  convex_to_fill_mesh = cv_mesh_ptr.release();
  convex_to_fill_geometry = cv_geometry_ptr.release();
  preprocess_mesh(convex_to_fill_mesh, convex_to_fill_geometry, triangulate || std::strcmp(poly_str.c_str(), "gomboc") == 0);
  convex_to_fill_geometry->inputVertexPositions *= 3.;
  // convex_to_fill_geometry->inputVertexPositions += Vector3::constant(0.5);

  // Register the mesh with polyscope
  polyscope::SurfaceMesh *psHullFillMesh = polyscope::registerSurfaceMesh(
    "fillable hull",
    convex_to_fill_geometry->inputVertexPositions, convex_to_fill_mesh->getFaceVertexList(),
    polyscopePermutations(*convex_to_fill_mesh));
  psHullFillMesh->setTransparency(0.75);
}


void generate_polyhedron_example(std::string poly_str, bool triangulate = false){
    // readManifoldSurfaceMesh()
  std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(poly_str);
  mesh = mesh_ptr.release();
  geometry = geometry_ptr.release();
  preprocess_mesh(mesh, geometry, triangulate || std::strcmp(poly_str.c_str(), "gomboc") == 0);
  first_time = true;
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
  // boundary_builder->print_area_of_boundary_loops();
  // printf("precomputes took %d\n", clock() - t1);
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
  // printf("GM visuals took %d\n", clock() - t1);
  t1 = clock();
  // initialize_boundary_builder();
  // printf("boundary took %d\n", clock() - t1);
  t1 = clock();
  if (draw_boundary_patches)
    draw_stable_patches_on_gauss_map(gm_is_drawn);
  // printf("boundary visuals took %d\n", clock() - t1);
  // initiate_markov_model();
  // visualize_sudo_faces();
}

void update_hull_of_hull(VertexData<Vector3> positions){
  std::vector<std::vector<size_t>> hull_faces; 
  std::vector<size_t> hull_vertex_mapping;
  std::vector<Vector3> hull_poses; // redundant, but helps with keeping this function clean
  std::tie(hull_faces, hull_vertex_mapping, hull_poses) = get_convex_hull(positions);

  forwardSolver->hullMesh = new ManifoldSurfaceMesh(hull_faces);
  forwardSolver->hullGeometry = new VertexPositionGeometry(*forwardSolver->hullMesh);
  for (Vertex v: forwardSolver->hullMesh->vertices())
    forwardSolver->hullGeometry->inputVertexPositions[v] = hull_poses[v.getIndex()];
}

// hull update stuff
double hull_update_line_search(VertexData<Vector3> grad){
  printf(" ---- hull update line search ----\n");
  double grad_norm2 = 0.;
  for (Vector3 v3: grad.toVector())
    grad_norm2 += v3.norm2();
  
  VertexData<Vector3> initial_poses = forwardSolver->hullGeometry->inputVertexPositions;
  std::vector<std::vector<size_t>> old_face_list = forwardSolver->hullMesh->getFaceVertexList();
  ManifoldSurfaceMesh *tmp_mesh = new ManifoldSurfaceMesh(forwardSolver->hullMesh->getFaceVertexList());
  //TODO: revert the current solver
  VertexData<Vector3> tmp_poses = forwardSolver->hullGeometry->inputVertexPositions;
  Forward3DSolver *tmp_solver = new Forward3DSolver(tmp_mesh, forwardSolver->hullGeometry,
                                                    forwardSolver->get_G());
  BoundaryBuilder *tmp_builder = new BoundaryBuilder(tmp_solver);

  tmp_solver->update_convex_hull(false); // changes tmp folver's hull mesh ..
  tmp_solver->updated = false;
  tmp_solver->initialize_pre_computes();
  tmp_builder->build_boundary_normals();
  double s_min_dice_energy = tmp_builder->get_fair_dice_energy(fair_sides_count);
  printf(" current fair dice energy %f\n", s_min_dice_energy);
  double s_max = step_size3,
         s_min = 0.,
         decay = 0.8,
         _armijo_const = 1e-4;
  int max_iters = 100;
  double s = s_max/2.; //

  double tmp_fair_dice_energy;
  int j;
  for (j = 0; j < max_iters; j++) {
      // update stuff
      // printf(" ^^ at line search iter: %d  s = %f\n", j, s);
      tmp_solver->inputGeometry->inputVertexPositions = initial_poses + s * grad; // old solver's hull mesh
      tmp_solver->update_convex_hull(false); // changes tmp folver's "hull" mesh ..
      tmp_solver->updated = false;
      tmp_solver->initialize_pre_computes();
      tmp_builder->build_boundary_normals();
      
      tmp_fair_dice_energy = tmp_builder->get_fair_dice_energy(fair_sides_count);
      // printf("  *** temp fair dice energy %d: %f\n", j, tmp_fair_dice_energy);
      if (tmp_fair_dice_energy <= s_min_dice_energy - _armijo_const * s * grad_norm2){
          break; //  x new is good
      }
      else
          s *= decay;
  }
  printf("line search for dice ended at iter %d, s: %f, \n \t\t\t\t\t fnew: %f \n", j, s, tmp_fair_dice_energy);
  forwardSolver->hullGeometry->inputVertexPositions = initial_poses; // back to normal
  return s;
} 


void take_uni_mass_opt_vertices_step(){
  printf("finding uniform mass\n");
  forwardSolver->set_uniform_G();
  printf("initalize precomputes\n");
  forwardSolver->initialize_pre_computes();
  printf("finding derivatives\n");
  inverseSolver->find_uni_mass_d_pf_dv();
  printf(" uni mass derivatives found!\n");
  if (structured_opt){ // if updating the flow structure
    boundary_builder->build_boundary_normals(); // face-last-face is called 
    // boundary_builder->print_area_of_boundary_loops();
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
  double opt_step_size = hull_update_line_search(total_uni_mass_vertex_grads);

  polyscope::registerSurfaceMesh("pre-step hull", forwardSolver->hullGeometry->inputVertexPositions,
                                                  forwardSolver->hullMesh->getFaceVertexList());
  polyscope::registerSurfaceMesh("pre-step mesh", forwardSolver->inputGeometry->inputVertexPositions,
                                                  forwardSolver->inputMesh->getFaceVertexList())->setEnabled(false);
  
  // PC
  // polyscope::PointCloud* pre_pc = polyscope::registerPointCloud("pre-step hull points", forwardSolver->hullGeometry->inputVertexPositions);
  // pre_pc->setPointRadius(vis_utils.gm_pt_radi * 1., false);
  // pre_pc->setPointColor(glm::vec3({0.8,0.1, 0.1}));
  // pre_pc->setEnabled(false);
  
  if (deform_after){
    if (true){ // new bending stuff
      // forwardSolver's hull mesh and geometry updated
      update_hull_of_hull(forwardSolver->hullGeometry->inputVertexPositions 
                          + total_uni_mass_vertex_grads * opt_step_size); // updates hull mesh and geometry
      polyscope::registerSurfaceMesh("pree-deform hull mesh", forwardSolver->hullGeometry->inputVertexPositions,
                                                              forwardSolver->hullMesh->getFaceVertexList())->setEnabled(false);

      // deforming
      deformationSolver = new DeformationSolver(mesh, inverseSolver->initial_geometry, // TODO: or use the ex inputGeometry?
                                                forwardSolver->hullMesh, forwardSolver->hullGeometry);   
      deformationSolver->CP_lambda = pow(10, CP_lambda_exp);
      deformationSolver->CP_mu = CP_mu;
      deformationSolver->filling_max_iter = filling_max_iter;
      DenseMatrix<double> new_points = deformationSolver->solve_for_bending(1);
      for (Vertex v: forwardSolver->inputMesh->vertices())
        forwardSolver->inputGeometry->inputVertexPositions[v] = vec_to_GC_vec3(new_points.row(v.getIndex()));
    }
    else if (false){ // old ARAP stuff
      printf("update hard correspondence\n");
      forwardSolver->update_hull_points_correspondence(forwardSolver->hullGeometry->inputVertexPositions 
                                                        + total_uni_mass_vertex_grads * opt_step_size, 
                                                        inverseSolver->initial_geometry->inputVertexPositions);
      
      VertexData<Vector3> trivial_updates;
      trivial_updates = inverseSolver->trivial_update_positions(total_uni_mass_vertex_grads * opt_step_size);
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
    }
  }
  else {
    VertexData<Vector3> trivial_updates;
    trivial_updates = inverseSolver->trivial_update_positions(total_uni_mass_vertex_grads * opt_step_size);
    forwardSolver->inputGeometry->inputVertexPositions = forwardSolver->inputGeometry->inputVertexPositions + trivial_updates;
  }
  polyscope::registerSurfaceMesh("post-deform pre-hull-update hull mesh", forwardSolver->hullGeometry->inputVertexPositions,
                                              forwardSolver->hullMesh->getFaceVertexList())->setEnabled(false);
  forwardSolver->update_convex_hull(with_hull_projection);
  // update hull with new positions
  
  polyscope::registerSurfaceMesh("post-deformation mesh", forwardSolver->inputGeometry->inputVertexPositions,
                                               forwardSolver->inputMesh->getFaceVertexList());
  polyscope::registerSurfaceMesh("post-deform & hull update hull mesh", forwardSolver->hullGeometry->inputVertexPositions,
                                              forwardSolver->hullMesh->getFaceVertexList())->setEnabled(false);
  
  // PC
  // polyscope::PointCloud* post_pc = polyscope::registerPointCloud("post-step hull points", forwardSolver->hullGeometry->inputVertexPositions);
  // pre_pc->setPointRadius(vis_utils.gm_pt_radi * 1., false);
  // pre_pc->setPointColor(glm::vec3({0.1,0.1, 0.8}));
  forwardSolver->set_uniform_G();
  G = forwardSolver->get_G();
  update_solver_and_boundaries();
  boundary_builder->print_area_of_boundary_loops();
  printf(" fair dice energy: %f\n", boundary_builder->get_fair_dice_energy(fair_sides_count));
  // update_visuals_with_G();
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
              init_visuals();
              // recolor_faces = true;
              // color_faces();

              // //
              // visualize_gauss_map();//
              // G = find_center_of_mass(*forwardSolver->inputMesh, *forwardSolver->inputGeometry).first;
              // update_solver_and_boundaries();
              // update_visuals_with_G();
          }
          if (is_selected)
              ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
      }
      ImGui::EndCombo();
  }

  if (ImGui::BeginCombo("##combo2", all_polygons_current_item2.c_str())){ // The second parameter is the label previewed before opening the combo.
      for (std::string tmp_str: all_polyhedra_items){ // This enables not having to have a const char* arr[]. Or maybe I'm just a noob.
          bool is_selected = (all_polygons_current_item2 == tmp_str.c_str()); // You can store your selection however you want, outside or inside your objects
          if (ImGui::Selectable(tmp_str.c_str(), is_selected)){ // selected smth
              all_polygons_current_item2 = tmp_str;
              init_convex_shape_to_fill(all_polygons_current_item2);
              deformationSolver = new DeformationSolver(mesh, geometry, convex_to_fill_mesh, convex_to_fill_geometry);
              // 
          }
          if (is_selected)
              ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
      }
      ImGui::EndCombo();
  }

  if (ImGui::Button("deform into convex shape")){
    init_convex_shape_to_fill(all_polygons_current_item2);
    deformationSolver = new DeformationSolver(mesh, geometry, convex_to_fill_mesh, convex_to_fill_geometry);   
    deformationSolver->CP_lambda = pow(10, CP_lambda_exp);
    deformationSolver->filling_max_iter = filling_max_iter;
    deformationSolver->CP_mu = CP_mu;
    deformationSolver->one_time_CP_assignment = one_time_CP_assignment;
    deformationSolver->solve_for_bending(1);
  }
  if (ImGui::Checkbox("one time CP assignment", &one_time_CP_assignment)) deformationSolver->one_time_CP_assignment = one_time_CP_assignment;
  if (ImGui::SliderFloat("CP lambda log10/initial value", &CP_lambda_exp, 0., 10.)) deformationSolver->CP_lambda = pow(10, CP_lambda_exp);
  if (ImGui::SliderFloat("CP mu (growth per iter)", &CP_mu, 1., 1.5)) deformationSolver->CP_mu = CP_mu;
  if (ImGui::SliderInt("filling iters", &filling_max_iter, 0., 200.0)) deformationSolver->filling_max_iter = filling_max_iter;

  if (ImGui::Button("uniform mass G")){
    G = find_center_of_mass(*forwardSolver->inputMesh, *forwardSolver->inputGeometry).first;
    update_solver_and_boundaries();
    update_visuals_with_G();
  }
  if (ImGui::Checkbox("draw artificial R3 boundaries", &test_guess)) draw_stable_patches_on_gauss_map();
  
  if (ImGui::Button("take fair step (joint gradient)")) {
    take_uni_mass_opt_vertices_step();
  }
  if (ImGui::SliderFloat("starting step size", &step_size3, 0., 1.0));
  if (ImGui::SliderInt("fair dice side count", &fair_sides_count, 1, 10));

  if (ImGui::Checkbox("deform after", &deform_after));
  if (ImGui::Checkbox("line search for probability gradients", &line_search));
  
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
  init_visuals();

  // init_convex_shape_to_fill(all_polygons_current_item2);
  // deformationSolver = new DeformationSolver(mesh, geometry, convex_to_fill_mesh, convex_to_fill_geometry);
              
  // convex_hull(forwardSolver->hullGeometry->inputVertexPositions);
  // build the solver

  // Initialize polyscope
  polyscope::init();
  // Set the callback function
  polyscope::state::userCallback = myCallback;
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  polyscope::view::upDir = polyscope::view::UpDir::NegZUp;
  
  // Give control to the polyscope gui
  polyscope::show();
  

  return EXIT_SUCCESS;
}
