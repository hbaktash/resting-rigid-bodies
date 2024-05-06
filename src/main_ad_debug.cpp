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
#include "inv_design.h"
// #include "igl/arap.h"
// #include "optimization.h"

#define ANSI_FG_MAGENTA "\x1b[35m"
#define ANSI_FG_YELLOW "\x1b[33m"
#define ANSI_FG_GREEN "\x1b[32m"
#define ANSI_FG_WHITE "\x1b[37m"
#define ANSI_FG_RED "\x1b[31m"
#define ANSI_RESET "\x1b[0m"


using namespace geometrycentral;
using namespace geometrycentral::surface;


// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr, cv_mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr, cv_geometry_ptr;
ManifoldSurfaceMesh *mesh;
VertexPositionGeometry *geometry, *deformed_geometry;
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
     first_time = false,
     curvature_weighted_CP = false;
bool do_remesh = true;
float scale_for_remesh = 1.003;
polyscope::PointCloud *test_pc;

int fair_sides_count = 6; // for optimization
bool do_sobolev_dice_grads = true,
     use_autodiff_for_dice_grad = true;
float sobolev_lambda = 2.,
      sobolev_lambda_decay = 0.95;
int sobolev_p = 2;
// optimization stuff
InverseSolver* inverseSolver;
float step_size = 0.01,
      step_size2 = 0.01,
      step_size3 = 0.1; // 0.167 small_bunny
float stable_normal_update_thresh = -1;

// deformation
bool v2_dice_animate = false;
float dice_search_decay = 0.95;
int hull_opt_steps = 20;

// static double xs2[20], ys2[20];

// test energies stuff
float row0[3] = {2.,0.,0.}; 
float row1[3] = {0.,2.,0.}; 
float row2[3] = {0.,0.,2.}; 

// example choice
std::vector<std::string> all_polyhedra_items = {std::string("tet"), std::string("tet2"), std::string("cube"), std::string("tilted cube"), std::string("sliced tet"), std::string("Conway spiral 4"), std::string("oloid"), std::string("fox"), std::string("small_bunny"), std::string("bunnylp"), std::string("kitten"), std::string("double-torus"), std::string("soccerball"), std::string("cowhead"), std::string("bunny"), std::string("gomboc"), std::string("mark_gomboc")};
std::string all_polygons_current_item = "small_bunny",
            all_polygons_current_item2 = "tet";
static const char* all_polygons_current_item_c_str = "bunnylp";


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
  preprocess_mesh(mesh, geometry, triangulate || std::strcmp(poly_str.c_str(), "gomboc") == 0, do_remesh, scale_for_remesh);

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
  // boundary_builder->build_boundary_normals();
  printf("  \n ************************************************ no AD \n");
  boundary_builder->build_boundary_normals();
  printf("  \n ************************************************  AD  \n");
  // autodiff::MatrixX3var poses_ad = vertex_data_to_matrix(forwardSolver->hullGeometry->inputVertexPositions);
  // autodiff::Vector3var G_ad = vec32vec(forwardSolver->get_G());
  boundary_builder = new BoundaryBuilder(forwardSolver);
  boundary_builder->build_boundary_normals_for_autodiff(true); // ( poses_ad, G_ad, ) autodiff; generate_gradients = true
  // for (Face f: forwardSolver->hullMesh->faces()){
  //   if (forwardSolver->face_last_face[f] == f){
  //     polyscope::getSurfaceMesh("init hull mesh")->addVertexVectorQuantity("ad grads", boundary_builder->df_dv_grads_ad[f])->setEnabled(true);
  //     break;
  //   }
  // }
}


void update_visuals_with_G(Forward3DSolver *tmp_solver = nullptr, BoundaryBuilder *bnd_builder = boundary_builder){
  if (tmp_solver != nullptr){
    tmp_solver->updated = false;
    tmp_solver->initialize_pre_computes();
    vis_utils.forwardSolver = tmp_solver;
  }
  
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
double hull_update_line_search(VertexData<Vector3> grad, Forward3DSolver fwd_solver, bool frozen_G){
  printf(" ---- hull update line search ----\n");
  // if (polyscope::hasSurfaceMesh("init hull mesh")) polyscope::getSurfaceMesh("init hull mesh")->setEnabled(true);
  double grad_norm2 = 0.;
  for (Vector3 v3: grad.toVector())
    grad_norm2 += v3.norm2();
  
  VertexData<Vector3> initial_poses = fwd_solver.hullGeometry->inputVertexPositions;
  std::vector<std::vector<size_t>> old_face_list = fwd_solver.hullMesh->getFaceVertexList();
  //TODO: revert the current solver
  ManifoldSurfaceMesh *tmp_mesh = new ManifoldSurfaceMesh(old_face_list);
  VertexPositionGeometry *tmp_geo = new VertexPositionGeometry(*tmp_mesh, vertex_data_to_matrix(initial_poses));
  Forward3DSolver *tmp_solver = new Forward3DSolver(tmp_mesh, tmp_geo,
                                                    fwd_solver.get_G(), false);
  tmp_solver->updated = false;
  tmp_solver->initialize_pre_computes();

  BoundaryBuilder *tmp_builder = new BoundaryBuilder(tmp_solver);
  tmp_builder->build_boundary_normals();

  double s_min_dice_energy = tmp_builder->get_fair_dice_energy(fair_sides_count);
  printf(" current fair dice energy: %f\n", s_min_dice_energy);
  double s_max = step_size3,
         s_min = 0.,
         _armijo_const = 0.; //1e-4;
  int max_iters = 400;
  double s = s_max; //

  double tmp_fair_dice_energy;
  int j;
  bool found_smth_optimal = false;
  for (j = 0; j < max_iters; j++) {
      // update stuff
      // printf(" ^^ at line search iter: %d  s = %f\n", j, s);
      auto [new_hull_mesh, new_hull_geo] = get_convex_hull_mesh(initial_poses + s * grad); // changes tmp folver's "hull" mesh ..
      tmp_solver = new Forward3DSolver(new_hull_mesh, new_hull_geo, tmp_solver->get_G(), false);
      if (!frozen_G)
        tmp_solver->set_uniform_G();
      tmp_solver->updated = false;
      tmp_solver->initialize_pre_computes();
      tmp_builder = new BoundaryBuilder(tmp_solver);
      tmp_builder->build_boundary_normals();
      
      tmp_fair_dice_energy = tmp_builder->get_fair_dice_energy(fair_sides_count);
      // printf("  *** temp fair dice energy %d: %f\n", j, tmp_fair_dice_energy);
      if (tmp_fair_dice_energy <= s_min_dice_energy - _armijo_const * s * grad_norm2){
        found_smth_optimal = true;
        break; //  x new is good
      }
      else
          s *= dice_search_decay;
  }
  s = found_smth_optimal ? s : 0.;
  printf("line search for dice ended at iter %d, s: %.10f, \n \t\t\t\t\t fnew: %f \n", j, s, tmp_fair_dice_energy);
  return s;
} 

void test_approx_vs_ad_grads(){
    // forwardSolver->set_uniform_G();
    ManifoldSurfaceMesh *hull_mesh = new ManifoldSurfaceMesh(forwardSolver->hullMesh->getFaceVertexList());
    VertexPositionGeometry *hull_geo = new VertexPositionGeometry(*hull_mesh, vertex_data_to_matrix(forwardSolver->hullGeometry->inputVertexPositions));
    Forward3DSolver *tmp_solver = new Forward3DSolver(hull_mesh, hull_geo, G, false);
    if (frozen_G){
    auto GV_pair = find_center_of_mass(*forwardSolver->inputMesh, *forwardSolver->inputGeometry);
    tmp_solver->set_G(GV_pair.first);
    tmp_solver->volume = GV_pair.second;
    }
    else {
    tmp_solver->set_uniform_G(); // frozen G matters here or not ???
    }
    printf("initalize precomputes\n");
    tmp_solver->initialize_pre_computes();

    if (polyscope::hasSurfaceMesh("fillable hull")) polyscope::getSurfaceMesh("fillable hull")->setEnabled(false);
    if (polyscope::hasSurfaceMesh("temp sol")) polyscope::getSurfaceMesh("temp sol")->setEnabled(false);
    if (polyscope::hasSurfaceMesh("init hull mesh")) polyscope::getSurfaceMesh("init hull mesh")->setEnabled(false);
    if (polyscope::hasSurfaceMesh("init input mesh")) polyscope::getSurfaceMesh("init input mesh")->setEnabled(false);
    BoundaryBuilder *tmp_bnd_builder = new BoundaryBuilder(tmp_solver);
    tmp_bnd_builder->build_boundary_normals_for_autodiff(true); // (poses_ad, G_ad, ) // autodiff; generate_gradients = true
    InverseSolver *tmp_inv_solver = new InverseSolver(tmp_bnd_builder);

    printf("finding vertex derivatives\n");
    printf("building normals bounds with AD\n");
    
    update_visuals_with_G(tmp_solver, tmp_bnd_builder);
    std::cout << ANSI_FG_YELLOW << " fair dice energy: " << tmp_bnd_builder->get_fair_dice_energy(fair_sides_count) << ANSI_RESET << "\n";

    auto pip2psmesh = polyscope::registerSurfaceMesh("pipe2 tmp sol", tmp_solver->inputGeometry->inputVertexPositions,
                                                  tmp_solver->inputMesh->getFaceVertexList());
    pip2psmesh->setEnabled(true);
    pip2psmesh->setTransparency(0.7);

    // dice energy policy
    if (first_time || always_update_structure){
      tmp_inv_solver->flow_structure = tmp_solver->face_last_face; // face last face is built in build_boundary_normals...
      if (first_time) stable_normal_update_thresh = -1.;
      else stable_normal_update_thresh = 0.1;
      first_time = false;
    }

    // fetching approx grad vs real grads
    printf("here\n");
    tmp_inv_solver->find_uni_mass_d_pf_dv(false, frozen_G);
    VertexData<Vector3> approx_dice_energy_grads = tmp_inv_solver->find_uni_mass_total_vertex_grads(fair_sides_count,
                                                                                               structured_opt, stable_normal_update_thresh);
    printf("here2\n");
    tmp_inv_solver->find_uni_mass_d_pf_dv(true, frozen_G);
    VertexData<Vector3> dice_energy_grads = tmp_inv_solver->find_uni_mass_total_vertex_grads(fair_sides_count,
                                                                                        structured_opt, stable_normal_update_thresh);
    printf("registering 3\n");
    
    // DEBUG
    // for (Vertex v: tmp_solver->hullMesh->vertices())
    //   std::cout << "approx grad: " << approx_dice_energy_grads[v] << "  \nreal grad: " << dice_energy_grads[v] << "\n";
    
    pip2psmesh->addVertexVectorQuantity("ad total grads", dice_energy_grads)->setEnabled(true);
    pip2psmesh->addVertexVectorQuantity("approx total grads", approx_dice_energy_grads)->setEnabled(true);
    printf("registered\n");
    // polyscope::frameTick();
}

void version2_dice_pipeline(size_t step_count = 1){
  // forwardSolver->set_uniform_G();
  ManifoldSurfaceMesh *hull_mesh = new ManifoldSurfaceMesh(forwardSolver->hullMesh->getFaceVertexList());
  VertexPositionGeometry *hull_geo = new VertexPositionGeometry(*hull_mesh, vertex_data_to_matrix(forwardSolver->hullGeometry->inputVertexPositions));
  Forward3DSolver *tmp_solver = new Forward3DSolver(hull_mesh, hull_geo, G, false);
  if (frozen_G){
    auto GV_pair = find_center_of_mass(*forwardSolver->inputMesh, *forwardSolver->inputGeometry);
    tmp_solver->set_G(GV_pair.first);
    tmp_solver->volume = GV_pair.second;
  }
  else {
    tmp_solver->set_uniform_G(); // frozen G matters here or not ???
  }
  printf("initalize precomputes\n");
  tmp_solver->initialize_pre_computes();
  
  if (polyscope::hasSurfaceMesh("fillable hull")) polyscope::getSurfaceMesh("fillable hull")->setEnabled(false);
  if (polyscope::hasSurfaceMesh("temp sol")) polyscope::getSurfaceMesh("temp sol")->setEnabled(false);
  if (polyscope::hasSurfaceMesh("init hull mesh")) polyscope::getSurfaceMesh("init hull mesh")->setEnabled(false);
  if (polyscope::hasSurfaceMesh("init input mesh")) polyscope::getSurfaceMesh("init input mesh")->setEnabled(false);
  auto pip2psmesh = polyscope::registerSurfaceMesh("pipe2 tmp sol", tmp_solver->inputGeometry->inputVertexPositions,
                                                  tmp_solver->inputMesh->getFaceVertexList());
  pip2psmesh->setEdgeWidth(2.);
  polyscope::screenshot(false);
  BoundaryBuilder *tmp_bnd_builder;
  InverseSolver *tmp_inv_solver;
  tmp_bnd_builder = new BoundaryBuilder(tmp_solver);
  tmp_inv_solver = new InverseSolver(tmp_bnd_builder);

  double current_sobolev_lambda = sobolev_lambda;
  for (size_t iter = 0; iter < step_count; iter++){

    // tmp_bnd_builder->build_boundary_normals();
    // update_visuals_with_G();
    // printf(" vertex derivatives found!\n");
    printf("finding vertex derivatives\n");
    if (!use_autodiff_for_dice_grad)
      tmp_bnd_builder->build_boundary_normals(); // face-last-face is called 
    else {
      // autodiff::MatrixX3var poses_ad = vertex_data_to_matrix(tmp_solver->hullGeometry->inputVertexPositions);
      // autodiff::Vector3var G_ad = vec32vec(tmp_solver->get_G());
      printf("building normals bounds with AD\n");
      tmp_bnd_builder->build_boundary_normals_for_autodiff(true); // (poses_ad, G_ad, ) // autodiff; generate_gradients = true
    }
    update_visuals_with_G(tmp_solver, tmp_bnd_builder);
    std::cout << ANSI_FG_YELLOW << " fair dice energy iter "<< iter<< "  f: " << tmp_bnd_builder->get_fair_dice_energy(fair_sides_count) << ANSI_RESET << "\n";

    auto pip2psmesh = polyscope::registerSurfaceMesh("pipe2 tmp sol", tmp_solver->inputGeometry->inputVertexPositions,
                                                  tmp_solver->inputMesh->getFaceVertexList());

    // boundary_builder->print_area_of_boundary_loops();
    if (first_time || always_update_structure){
      tmp_inv_solver->flow_structure = tmp_solver->face_last_face; // face last face is built in build_boundary_normals...
      if (first_time) stable_normal_update_thresh = -1.;
      else stable_normal_update_thresh = 0.1;
      first_time = false;
    }

    // checking approx grad vs real grads
    tmp_inv_solver->find_uni_mass_d_pf_dv(false, frozen_G);
    VertexData<Vector3> approx_dice_energy_grads = tmp_inv_solver->find_uni_mass_total_vertex_grads(fair_sides_count,
                                                                                               structured_opt, stable_normal_update_thresh);
    
    tmp_inv_solver->find_uni_mass_d_pf_dv(use_autodiff_for_dice_grad, frozen_G);
    VertexData<Vector3> dice_energy_grads = tmp_inv_solver->find_uni_mass_total_vertex_grads(fair_sides_count,
                                                                                        structured_opt, stable_normal_update_thresh);
    auto old_tmp_geo = tmp_solver->inputGeometry->inputVertexPositions;
    auto old_tmp_mesh = tmp_solver->inputMesh->getFaceVertexList();
    
    pip2psmesh->setEnabled(true);
    pip2psmesh->setTransparency(0.7);
    // DEBUG
    polyscope::getSurfaceMesh("pipe2 tmp sol")->addVertexVectorQuantity("ad total grads", dice_energy_grads)->setEnabled(true);
    polyscope::getSurfaceMesh("pipe2 tmp sol")->addVertexVectorQuantity("approx total grads", approx_dice_energy_grads)->setEnabled(true);
    polyscope::show();
    // for (Face f: tmp_solver->hullMesh->faces()){
    //   if (tmp_solver->face_last_face[f] == f){
    //     polyscope::getSurfaceMesh("pipe2 tmp sol")->addVertexVectorQuantity("ad grads", tmp_bnd_builder->df_dv_grads_ad[f])->setEnabled(true);
    //     polyscope::getSurfaceMesh("pipe2 tmp sol")->addVertexVectorQuantity("approx grads", approx_dpf_dv[f])->setEnabled(true);
    //     polyscope::show();
    //   }
    // }

    if (do_sobolev_dice_grads){
      
      current_sobolev_lambda *= sobolev_lambda_decay;
      dice_energy_grads = tmp_inv_solver->sobolev_diffuse_gradients(dice_energy_grads, current_sobolev_lambda, sobolev_p);
      
      // std::cout << "energies: " << energies[0][iter] << " " << energies[1][iter] << " " << energies[2][iter] << "\n";
      pip2psmesh->addVertexVectorQuantity(" sobolev diffused dice grads", dice_energy_grads)->setEnabled(true);
    }

    pip2psmesh->setSurfaceColor({0.1,0.9,0.1});
    pip2psmesh->setEdgeWidth(2.);
    polyscope::screenshot(false);
    polyscope::frameTick();
    // polyscope::show();

    double opt_step_size = hull_update_line_search(dice_energy_grads, *tmp_solver, frozen_G);
    std::cout << ANSI_FG_RED << "opt step size: " << opt_step_size << ANSI_RESET << "\n";
    auto [new_hull_mesh, new_hull_geo] = get_convex_hull_mesh(tmp_solver->hullGeometry->inputVertexPositions + opt_step_size * dice_energy_grads);
    tmp_solver = new Forward3DSolver(new_hull_mesh, new_hull_geo, tmp_solver->get_G(), false); 
    // TODO: is this ok?
  
    // visuals
    if (!frozen_G){
      tmp_solver->set_uniform_G();
    }
    tmp_solver->updated = false;
    tmp_solver->initialize_pre_computes();
    tmp_bnd_builder = new BoundaryBuilder(tmp_solver);
    tmp_inv_solver = new InverseSolver(tmp_bnd_builder);
  }
  tmp_bnd_builder->build_boundary_normals();
  update_visuals_with_G(tmp_solver, tmp_bnd_builder);
  tmp_bnd_builder->print_area_of_boundary_loops();

  //DEBUG
  if (polyscope::hasSurfaceMesh("init hull mesh")) polyscope::getSurfaceMesh("init hull mesh")->setEnabled(true);

  polyscope::removeAllStructures();
  polyscope::SurfaceMesh *psHullFillMesh = polyscope::registerSurfaceMesh(
    "fillable hull", tmp_solver->inputGeometry->inputVertexPositions, tmp_solver->inputMesh->getFaceVertexList());
  psHullFillMesh->setTransparency(0.45);
  psHullFillMesh->setEdgeWidth(2.);
  psHullFillMesh->setEnabled(true);
  gm_is_drawn = false;
  
  pre_deform_G = tmp_solver->get_G();
  std::cout << "pre-deform G:" << tmp_solver->get_G() << "\n";
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
  if (ImGui::Checkbox("do remeshing", &do_remesh));
  if (ImGui::SliderFloat("remesh edge len scale", &scale_for_remesh, 0.1, 4.));
  if (ImGui::BeginCombo("##combo1", all_polygons_current_item.c_str())){ // The second parameter is the label previewed before opening the combo.
      for (std::string tmp_str: all_polyhedra_items){ // This enables not having to have a const char* arr[]. Or maybe I'm just a noob.
          bool is_selected = (all_polygons_current_item == tmp_str.c_str()); // You can store your selection however you want, outside or inside your objects
          if (ImGui::Selectable(tmp_str.c_str(), is_selected)){ // selected smth
              all_polygons_current_item = tmp_str;
              generate_polyhedron_example(all_polygons_current_item);
              update_solver();
              init_visuals();

              // //
              visualize_gauss_map();//
              G = find_center_of_mass(*forwardSolver->inputMesh, *forwardSolver->inputGeometry).first;
              update_solver_and_boundaries();
              boundary_builder->print_area_of_boundary_loops();
              update_visuals_with_G(forwardSolver, boundary_builder);
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
  if (ImGui::Checkbox("draw artificial R3 boundaries", &test_guess)) {
    if (test_guess)
      draw_stable_patches_on_gauss_map();
    else{
      polyscope::getPointCloud("test point cloud")->setEnabled(false);
      polyscope::getCurveNetwork("stable regions on polyhedra")->setEnabled(false);
    }
  }

  
  if (ImGui::SliderFloat("starting step size", &step_size3, 0., 1.0));
  if (ImGui::SliderInt("fair dice side count", &fair_sides_count, 1, 10));
  if (ImGui::Checkbox("frozen G", &frozen_G));

  if (ImGui::SliderFloat("dice energy step decay", &dice_search_decay, 0., 1.));
  if (ImGui::SliderInt("hull optimize step count", &hull_opt_steps, 1, 200));
  if (ImGui::SliderFloat("stable normal update thresh", &stable_normal_update_thresh, 0., 4.0));
  if (ImGui::Checkbox("Sobolev pre-condition", &do_sobolev_dice_grads));
  if (do_sobolev_dice_grads){
    if (ImGui::SliderFloat("Sobolev lambda", &sobolev_lambda, 0., 5.));
    if (ImGui::SliderFloat("Sobolev lambda decay", &sobolev_lambda_decay, 0., 1));
    if (ImGui::SliderInt("Sobolev power", &sobolev_p, 1, 5));
  }
  if (ImGui::Checkbox("autodiff dice grads", &use_autodiff_for_dice_grad));
  if (ImGui::Button("test ad grads vs approx grads")) {
    test_approx_vs_ad_grads();
  }
//   if (ImGui::Button("optimize hull (w.r.t. dice energy)")) {
//     v2_dice_animate = true;
//     if (polyscope::hasSurfaceMesh("v2pipeline final mesh")) polyscope::removeSurfaceMesh("v2pipeline final mesh");
//     if (polyscope::hasSurfaceMesh("v2pipeline final hull")) polyscope::removeSurfaceMesh("v2pipeline final hull");
//   }
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

  
  // polyscope::show();
  while (true) {
    if (v2_dice_animate){
      // polyscope::removeAllStructures();
      version2_dice_pipeline(hull_opt_steps);
      v2_dice_animate = false;
    }
    polyscope::frameTick();
  }

  return EXIT_SUCCESS;
}
