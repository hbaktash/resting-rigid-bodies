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
// #include "optimization.h"
#include "implot.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr, cv_mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr, cv_geometry_ptr;
ManifoldSurfaceMesh* mesh, *convex_to_fill_mesh;
VertexPositionGeometry* geometry, *convex_to_fill_geometry;
bool load_hull_from_file = true;
Vector3 G, // center of Mass
        pre_deform_G({-10.,-10.,-10.}),
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
bool test_guess = false;
bool compute_global_G_effect = true,
     deform_after = true,
     frozen_G = false,
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
bool do_sobolev_dice_grads = true;
float sobolev_lambda = 0.1,
      sobolev_lambda_decay = 0.95;
int sobolev_p = 2;
// optimization stuff
InverseSolver* inverseSolver;
float step_size = 0.01,
      step_size2 = 0.01,
      step_size3 = 0.1; // 0.167 small_bunny
float stable_normal_update_thresh = -1;
int ARAP_max_iters = 10;

// deformation
DeformationSolver *deformationSolver;
bool animate = false,
     v2_dice_animate = false;
float dice_search_decay = 0.95;

float bending_lambda_exps[2] = {1., 1.},
      membrane_lambda_exps[2] = {3., 3.},
      CP_lambda_exps[2] = {1., 7.},
      barrier_lambda_exps[2] = {-4., -8.},
      reg_lambda_exp = -3.,
      internal_p = 0.91,
      refinement_CP_threshold = 0.05;
int filling_max_iter = 10;
int hull_opt_steps = 50;

static const int MAX_FILL_ITERS = 300,
                 ENERGY_COUNT = 5;
int current_ENERGY_COUNT = 3;
int current_fill_iter = -1;
static float xs[MAX_FILL_ITERS];
const char* energy_names[ENERGY_COUNT] = {"bending", "membrane", "CP", "barrier", "total"};
static float **energies = new float*[ENERGY_COUNT];
float max_energy = 0.;

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


void init_convex_shape_to_fill(std::string poly_str, bool triangulate = true){
    // readManifoldSurfaceMesh()
  std::tie(cv_mesh_ptr, cv_geometry_ptr) = generate_polyhedra(poly_str);
  convex_to_fill_mesh = cv_mesh_ptr.release();
  convex_to_fill_geometry = cv_geometry_ptr.release();
  preprocess_mesh(convex_to_fill_mesh, convex_to_fill_geometry, triangulate || std::strcmp(poly_str.c_str(), "gomboc") == 0);

  // Register the mesh with polyscope
  polyscope::SurfaceMesh *psHullFillMesh = polyscope::registerSurfaceMesh(
    "fillable hull",
    convex_to_fill_geometry->inputVertexPositions, convex_to_fill_mesh->getFaceVertexList(),
    polyscopePermutations(*convex_to_fill_mesh));
  psHullFillMesh->setTransparency(0.35);
}

void initialize_deformation_params(DeformationSolver *deformation_solver){
  deformationSolver->dynamic_remesh = dynamic_remesh;
  deformation_solver->refinement_CP_threshold = refinement_CP_threshold;
  deformation_solver->filling_max_iter = filling_max_iter;  
  
  deformationSolver->init_bending_lambda = pow(10, bending_lambda_exps[0]);
  deformationSolver->final_bending_lambda = pow(10, bending_lambda_exps[1]);
  deformationSolver->init_membrane_lambda = pow(10, membrane_lambda_exps[0]);
  deformationSolver->final_membrane_lambda = pow(10, membrane_lambda_exps[1]);
  deformationSolver->init_CP_lambda = pow(10, CP_lambda_exps[0]);
  deformationSolver->final_CP_lambda = pow(10, CP_lambda_exps[1]);
  if (!use_QP_solver){
    deformationSolver->init_barrier_lambda = pow(10, barrier_lambda_exps[0]);
    deformationSolver->final_barrier_lambda = pow(10, barrier_lambda_exps[1]);
  }
  else {
    deformationSolver->init_barrier_lambda = 0.;
    deformationSolver->final_barrier_lambda = 0.;
  }
  deformation_solver->internal_growth_p = internal_p;
  if (use_reg)
    deformation_solver->reg_lambda = pow(10, reg_lambda_exp);
  else 
    deformation_solver->reg_lambda = 0.;
  deformation_solver->curvature_weighted_CP = curvature_weighted_CP;
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
  boundary_builder->build_boundary_normals();
  // boundary_builder->print_area_of_boundary_loops();
  // printf("precomputes took %d\n", clock() - t1);
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
double hull_update_line_search(VertexData<Vector3> grad, Forward3DSolver *fwd_solver, bool frozen_G){
  printf(" ---- hull update line search ----\n");
  // if (polyscope::hasSurfaceMesh("init hull mesh")) polyscope::getSurfaceMesh("init hull mesh")->setEnabled(true);
  double grad_norm2 = 0.;
  for (Vector3 v3: grad.toVector())
    grad_norm2 += v3.norm2();
  
  VertexData<Vector3> initial_poses = fwd_solver->hullGeometry->inputVertexPositions;
  std::vector<std::vector<size_t>> old_face_list = fwd_solver->hullMesh->getFaceVertexList();
  //TODO: revert the current solver
  ManifoldSurfaceMesh *tmp_mesh = new ManifoldSurfaceMesh(old_face_list);
  VertexPositionGeometry *tmp_geo = new VertexPositionGeometry(*tmp_mesh, vertex_data_to_matrix(initial_poses));
  Forward3DSolver *tmp_solver = new Forward3DSolver(tmp_mesh, tmp_geo,
                                                    fwd_solver->get_G(), false);
  BoundaryBuilder *tmp_builder = new BoundaryBuilder(tmp_solver);
  tmp_solver->updated = false;
  tmp_solver->initialize_pre_computes();
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
      tmp_solver = new Forward3DSolver(new_hull_mesh, new_hull_geo, fwd_solver->get_G(), false);
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
  printf("line search for dice ended at iter %d, s: %.10f, \n \t\t\t\t\t fnew: %f \n", j, s, tmp_fair_dice_energy);
  return found_smth_optimal ? s : 0.;
} 


// void take_uni_mass_opt_vertices_step(bool frozen_G = false){
//   if (!frozen_G){
//     printf("finding uniform mass\n");
//     forwardSolver->set_uniform_G();
//   }
//   printf("initalize precomputes\n");
//   forwardSolver->initialize_pre_computes();
//   printf("finding vertex derivatives\n");
//   inverseSolver->find_uni_mass_d_pf_dv(frozen_G);
//   printf(" vertex derivatives found!\n");
//   if (structured_opt){ // if updating the flow structure
//     boundary_builder->build_boundary_normals(); // face-last-face is called 
//     // boundary_builder->print_area_of_boundary_loops();
//     if (first_time || always_update_structure){
//       inverseSolver->flow_structure = forwardSolver->face_last_face;
//       if (first_time) stable_normal_update_thresh = -1.;
//       else stable_normal_update_thresh = 0.1;
//       first_time = false;
//     }
//   }
//   printf(" finding total per vertex derivatives\n");
//   VertexData<Vector3> total_uni_mass_vertex_grads = inverseSolver->find_uni_mass_total_vertex_grads(structured_opt, stable_normal_update_thresh);
//   printf(" total per vertex derivatives found!\n");
//   double opt_step_size = hull_update_line_search(total_uni_mass_vertex_grads, forwardSolver, frozen_G);

//   polyscope::registerSurfaceMesh("pre-step hull", forwardSolver->hullGeometry->inputVertexPositions,
//                                                   forwardSolver->hullMesh->getFaceVertexList())->setEnabled(false);
//   polyscope::registerSurfaceMesh("pre-step mesh", forwardSolver->inputGeometry->inputVertexPositions,
//                                                   forwardSolver->inputMesh->getFaceVertexList())->setEnabled(false);
  
//   // PC
//   // polyscope::PointCloud* pre_pc = polyscope::registerPointCloud("pre-step hull points", forwardSolver->hullGeometry->inputVertexPositions);
//   // pre_pc->setPointRadius(vis_utils.gm_pt_radi * 1., false);
//   // pre_pc->setPointColor(glm::vec3({0.8,0.1, 0.1}));
//   // pre_pc->setEnabled(false);
  
//   if (deform_after){
//     if (true){ // new bending stuff
//       // forwardSolver's hull mesh and geometry updated
//       update_hull_of_hull(forwardSolver->hullGeometry->inputVertexPositions 
//                           + total_uni_mass_vertex_grads * opt_step_size); // updates hull mesh and geometry
//       polyscope::registerSurfaceMesh("pree-deform hull mesh", forwardSolver->hullGeometry->inputVertexPositions,
//                                                               forwardSolver->hullMesh->getFaceVertexList())->setEnabled(false);

//       // deforming
//       deformationSolver = new DeformationSolver(mesh, inverseSolver->initial_geometry, // TODO: or use the ex inputGeometry?
//                                                 forwardSolver->hullMesh, forwardSolver->hullGeometry);   
//       initialize_deformation_params(deformationSolver);

//       DenseMatrix<double> new_points = deformationSolver->solve_for_bending(1);
//       for (Vertex v: forwardSolver->inputMesh->vertices())
//         forwardSolver->inputGeometry->inputVertexPositions[v] = vec_to_GC_vec3(new_points.row(v.getIndex()));
//     }
//     else if (false){ // old ARAP stuff
//       printf("update hard correspondence\n");
//       forwardSolver->update_hull_points_correspondence(forwardSolver->hullGeometry->inputVertexPositions 
//                                                         + total_uni_mass_vertex_grads * opt_step_size, 
//                                                         inverseSolver->initial_geometry->inputVertexPositions);
      
//       VertexData<Vector3> trivial_updates;
//       trivial_updates = inverseSolver->trivial_update_positions(total_uni_mass_vertex_grads * opt_step_size);
//       auto new_positions = forwardSolver->inputGeometry->inputVertexPositions + trivial_updates;
//       polyscope::registerSurfaceMesh("pre-ARAP mesh", new_positions,
//                                                   forwardSolver->inputMesh->getFaceVertexList())->setEnabled(false);
//       // update indexing for interior vertices


//       printf("updating interior positions\n");
//       VertexData<Vector3> updates;
//       if (forwardSolver->hullMesh->nVertices() == forwardSolver->inputMesh->nVertices()) // no interior vertices
//         updates = trivial_updates;
//       else {
//         // updates = inverseSolver->laplace_update_positions(new_positions);
//         updates = inverseSolver->ARAP_update_positions(new_positions);
//         // updates = inverseSolver->greedy_update_positions(total_uni_mass_vertex_grads * step_size3);
//         // updates = inverseSolver->diffusive_update_positions(total_uni_mass_vertex_grads * step_size3);
//       }
//       forwardSolver->inputGeometry->inputVertexPositions += updates;
//     }
//   }
//   else {
//     VertexData<Vector3> trivial_updates;
//     trivial_updates = inverseSolver->trivial_update_positions(total_uni_mass_vertex_grads * opt_step_size);
//     forwardSolver->inputGeometry->inputVertexPositions = forwardSolver->inputGeometry->inputVertexPositions + trivial_updates;
//   }
//   polyscope::registerSurfaceMesh("post-deform pre-hull-update hull mesh", forwardSolver->hullGeometry->inputVertexPositions,
//                                  forwardSolver->hullMesh->getFaceVertexList())->setEnabled(false);
//   forwardSolver->update_convex_hull(with_hull_projection);
//   // update hull with new positions
  
//   polyscope::registerSurfaceMesh("post-deformation mesh", forwardSolver->inputGeometry->inputVertexPositions,
//                                  forwardSolver->inputMesh->getFaceVertexList());
//   polyscope::registerSurfaceMesh("post-deform & hull update hull mesh", forwardSolver->hullGeometry->inputVertexPositions,
//                                  forwardSolver->hullMesh->getFaceVertexList())->setEnabled(false);
  
//   // PC
//   // polyscope::PointCloud* post_pc = polyscope::registerPointCloud("post-step hull points", forwardSolver->hullGeometry->inputVertexPositions);
//   // pre_pc->setPointRadius(vis_utils.gm_pt_radi * 1., false);
//   // pre_pc->setPointColor(glm::vec3({0.1,0.1, 0.8}));
//   if (!frozen_G){
//     forwardSolver->set_uniform_G();
//     G = forwardSolver->get_G();
//   }
//   update_solver_and_boundaries();
//   update_visuals_with_G();
//   boundary_builder->print_area_of_boundary_loops();
//   printf(" fair dice energy: %f\n", boundary_builder->get_fair_dice_energy(fair_sides_count));
//   // update_visuals_with_G();
// }


void animate_convex_fill_deformation(ManifoldSurfaceMesh *_mesh, VertexPositionGeometry *_old_geometry,
                                     ManifoldSurfaceMesh *_convex_mesh, VertexPositionGeometry *_convex_geometry){
  if (polyscope::hasSurfaceMesh("fillable hull")) polyscope::getSurfaceMesh("fillable hull")->setEnabled(true);
  deformationSolver = new DeformationSolver(_mesh, _old_geometry, _convex_mesh, _convex_geometry);   
  initialize_deformation_params(deformationSolver);
  animate = true;
}


void version2_dice_pipeline(size_t step_count = 1){
  // forwardSolver->set_uniform_G();
  ManifoldSurfaceMesh *hull_mesh = new ManifoldSurfaceMesh(forwardSolver->hullMesh->getFaceVertexList());
  VertexPositionGeometry *hull_geo = new VertexPositionGeometry(*hull_mesh, vertex_data_to_matrix(forwardSolver->hullGeometry->inputVertexPositions));
  Forward3DSolver *tmp_solver = new Forward3DSolver(hull_mesh, hull_geo, G, false);
  tmp_solver->set_uniform_G(); // frozen G matters here or not ???
  printf("initalize precomputes\n");
  tmp_solver->initialize_pre_computes();
  printf("finding vertex derivatives\n");
  
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

  max_energy = 0.;
  // cur = 2
  current_fill_iter = -1;
  current_ENERGY_COUNT = 3;
  energy_names[0] = "Sblv grad norm";
  energy_names[1] = "raw grad norm";
  energy_names[2] = "Sblv/raw grad norm";
  double current_sobolev_lambda = sobolev_lambda;
  for (size_t iter = 0; iter < step_count; iter++){

    tmp_inv_solver->find_uni_mass_d_pf_dv(frozen_G);
    // printf(" vertex derivatives found!\n");
    if (structured_opt){ // if updating the flow structure
      tmp_bnd_builder->build_boundary_normals(); // face-last-face is called 
      // boundary_builder->print_area_of_boundary_loops();
      if (first_time || always_update_structure){
        tmp_inv_solver->flow_structure = tmp_solver->face_last_face;
        if (first_time) stable_normal_update_thresh = -1.;
        else stable_normal_update_thresh = 0.1;
        first_time = false;
      }
    }
    VertexData<Vector3> vertex_grads = tmp_inv_solver->find_uni_mass_total_vertex_grads(fair_sides_count,
                                                                                        structured_opt, stable_normal_update_thresh);
    auto pip2psmesh = polyscope::registerSurfaceMesh("pipe2 tmp sol", tmp_solver->inputGeometry->inputVertexPositions,
                                                  tmp_solver->inputMesh->getFaceVertexList());
    auto old_tmp_geo = tmp_solver->inputGeometry->inputVertexPositions;
    auto old_tmp_mesh = tmp_solver->inputMesh->getFaceVertexList();
    
    pip2psmesh->setEnabled(true);
    pip2psmesh->setTransparency(0.7);
    pip2psmesh->addVertexVectorQuantity(" dice grads raw", vertex_grads)->setVectorLengthScale(0.04)->setEnabled(true);
    
    if (do_sobolev_dice_grads){
      
      current_fill_iter = iter;
      energies[1][iter] = tinyAD_flatten(vertex_data_to_matrix(vertex_grads)).norm();
      current_sobolev_lambda *= sobolev_lambda_decay;
      vertex_grads = tmp_inv_solver->sobolev_diffuse_gradients(vertex_grads, current_sobolev_lambda, sobolev_p);
      
      energies[0][iter] = tinyAD_flatten(vertex_data_to_matrix(vertex_grads)).norm();
      energies[2][iter] = energies[0][iter] / energies[1][iter];
      pip2psmesh->addVertexVectorQuantity(" sobolev diffused dice grads", vertex_grads)->setEnabled(true);
    }

    pip2psmesh->setSurfaceColor({0.1,0.9,0.1});
    pip2psmesh->setEdgeWidth(2.);
    polyscope::screenshot(false);
    polyscope::frameTick();
    // polyscope::show();

    double opt_step_size = hull_update_line_search(vertex_grads, tmp_solver, frozen_G);
    auto [new_hull_mesh, new_hull_geo] = get_convex_hull_mesh(tmp_solver->hullGeometry->inputVertexPositions + opt_step_size * vertex_grads);
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
  
    tmp_bnd_builder->build_boundary_normals();
    update_visuals_with_G(tmp_solver, tmp_bnd_builder);
    printf(" fair dice energy iter %d: %f\n", iter, tmp_bnd_builder->get_fair_dice_energy(fair_sides_count));
    // update_visuals_with_G();
  }
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
  convex_to_fill_mesh = tmp_solver->inputMesh;
  convex_to_fill_geometry = tmp_solver->inputGeometry;
  load_hull_from_file = false;
  
  pre_deform_G = tmp_solver->get_G();
  std::cout << "pre-deform G:" << tmp_solver->get_G() << "\n";
  std::cout << "*** Test ***\n";
  Forward3DSolver *tmp_solver2 = new Forward3DSolver(convex_to_fill_mesh, convex_to_fill_geometry, pre_deform_G, false);
  tmp_solver2->set_uniform_G();
  tmp_solver2->updated = false;
  tmp_solver2->initialize_pre_computes();
  BoundaryBuilder* tmp_bnd_builder2 = new BoundaryBuilder(tmp_solver2);
  tmp_bnd_builder2->build_boundary_normals();
  tmp_bnd_builder2->print_area_of_boundary_loops();
  if (draw_boundary_patches)
    draw_stable_patches_on_gauss_map(gm_is_drawn, tmp_bnd_builder2);
  std::cout << "dice energy check:" << tmp_bnd_builder2->get_fair_dice_energy(fair_sides_count) << "\n";
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
              // recolor_faces = true;
              // color_faces();

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
  if (ImGui::BeginCombo("##combo2", all_polygons_current_item2.c_str())){ // The second parameter is the label previewed before opening the combo.
      for (std::string tmp_str: all_polyhedra_items){ // This enables not having to have a const char* arr[]. Or maybe I'm just a noob.
          bool is_selected = (all_polygons_current_item2 == tmp_str.c_str()); // You can store your selection however you want, outside or inside your objects
          if (ImGui::Selectable(tmp_str.c_str(), is_selected)){ // selected smth
              all_polygons_current_item2 = tmp_str;
              init_convex_shape_to_fill(all_polygons_current_item2);
              deformationSolver = new DeformationSolver(mesh, geometry, convex_to_fill_mesh, convex_to_fill_geometry);
              load_hull_from_file = true;
              pre_deform_G = Vector3({-10.,-10.,-10.}); //
              // 
          }
          if (is_selected)
              ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
      }
      ImGui::EndCombo();
  }

  ImPlot::CreateContext();
  if (ImPlot::BeginPlot("CP + elastic energies", "iter", "energy")) {
      for (int i = 0; i < current_ENERGY_COUNT; i++)
        if (max_energy < energies[i][current_fill_iter]) max_energy = energies[i][current_fill_iter];
      ImPlot::SetupAxes("iter","energy");
      ImPlot::SetupAxisLimits(ImAxis_X1,0, current_fill_iter);
      ImPlot::SetupAxisLimits(ImAxis_Y1,0, max_energy);
      for (int i = 0; i < current_ENERGY_COUNT; i++)
        ImPlot::PlotLine(energy_names[i], xs, energies[i], current_fill_iter);
      ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle);
      ImPlot::EndPlot();
  }
  ImPlot::DestroyContext();

  if (ImGui::Button("deform into convex shape")){
    polyscope::removeAllStructures();
    if (load_hull_from_file)
      init_convex_shape_to_fill(all_polygons_current_item2);
    else {
      polyscope::SurfaceMesh *psHullFillMesh = polyscope::registerSurfaceMesh(
      "fillable hull",
      convex_to_fill_geometry->inputVertexPositions, convex_to_fill_mesh->getFaceVertexList(),
      polyscopePermutations(*convex_to_fill_mesh));
      psHullFillMesh->setTransparency(0.35);
    }
    animate_convex_fill_deformation(mesh, geometry, convex_to_fill_mesh, convex_to_fill_geometry);
    // TODO: check what should be done here to avoid the ugly bool trick
    // deformationSolver->solve_for_bending(1);
  }
  if (ImGui::SliderInt("filling iters", &filling_max_iter, 1, 300)) {
    deformationSolver->filling_max_iter = filling_max_iter;
  }
  if (ImGui::Checkbox("dynamic remesh", &dynamic_remesh));
  if (ImGui::SliderFloat("growth p", &internal_p, 0., 1.));
  if (ImGui::SliderFloat("refinement CP threshold ", &refinement_CP_threshold, 0., 1.));
  if (ImGui::InputFloat2("init/final bending log ", bending_lambda_exps)
    || ImGui::InputFloat2("init/final membrane log ", membrane_lambda_exps)
    || ImGui::InputFloat2("init/final CP log ", CP_lambda_exps)
    || ImGui::Checkbox("use QP solver ", &use_QP_solver));
  if(!use_QP_solver) 
    if (ImGui::InputFloat2("init/final barrier log ", barrier_lambda_exps));
  if(ImGui::Checkbox("curvature weighted", &curvature_weighted_CP));
  ImGui::Checkbox("use reg ", &use_reg);
  if(use_reg) 
    if (ImGui::SliderFloat("reg lambda; log10", &reg_lambda_exp, -5., 5.));
  

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

  // // V3 stuff 
  // if (ImGui::Button("take fair step (joint gradient)")) {
  //   take_uni_mass_opt_vertices_step(frozen_G);
  // }
  // if (ImGui::Checkbox("deform after", &deform_after));
  // if (ImGui::Checkbox("compute global G effect", &compute_global_G_effect)) inverseSolver->compute_global_G_effect = compute_global_G_effect;
  // if (ImGui::Checkbox("structured opt", &structured_opt));
  // if (ImGui::Checkbox("update structured at every step", &always_update_structure));
  // if (ImGui::Checkbox("decrease convex hull points", &with_hull_projection));

  if (ImGui::SliderFloat("dice energy step decay", &dice_search_decay, 0., 1.));
  if (ImGui::SliderInt("hull optimize step count", &hull_opt_steps, 1, 200));
  if (ImGui::SliderFloat("stable normal update thresh", &stable_normal_update_thresh, 0., 4.0));
  if (ImGui::Checkbox("Sobolev pre-condition", &do_sobolev_dice_grads));
  if (do_sobolev_dice_grads){
    if (ImGui::SliderFloat("Sobolev lambda", &sobolev_lambda, 0., 2.));
    if (ImGui::SliderFloat("Sobolev lambda decay", &sobolev_lambda_decay, 0., 1));
    if (ImGui::SliderInt("Sobolev power", &sobolev_p, 1, 5));
  }
  if (ImGui::Button("optimize hull (w.r.t. dice energy)")) {
    v2_dice_animate = true;
    if (polyscope::hasSurfaceMesh("v2pipeline final mesh")) polyscope::removeSurfaceMesh("v2pipeline final mesh");
    if (polyscope::hasSurfaceMesh("v2pipeline final hull")) polyscope::removeSurfaceMesh("v2pipeline final hull");
  }

  // Test stuff
  if (ImGui::InputFloat3("row0", row0)
    || ImGui::InputFloat3("row1", row1)
    || ImGui::InputFloat3("row2", row2));
  if (ImGui::Button("simple test")){
    DenseMatrix<double> A = DenseMatrix<double>::Zero(3,3);
    A(0,0) = row0[0]; A(0,1) = row0[1]; A(0,2) = row0[2];
    A(1,0) = row1[0]; A(1,1) = row1[1]; A(1,2) = row1[2];
    A(2,0) = row2[0]; A(2,1) = row2[1]; A(2,2) = row2[2];
    std::cout << "A matrix\n \t" << A << "\n";
    // auto tmp_def = new DeformationSolver(forwardSolver->inputMesh, forwardSolver->inputGeometry, convex_to_fill_mesh, convex_to_fill_geometry);   
    // initialize_deformation_params(tmp_def);
    // tmp_def->print_energies_after_transform(A);
    polyscope::removeAllStructures();
    convex_to_fill_mesh = forwardSolver->hullMesh;
    convex_to_fill_geometry = forwardSolver->hullGeometry;
    preprocess_mesh(convex_to_fill_mesh, convex_to_fill_geometry, true);

    // transform
    for (Vertex v: convex_to_fill_mesh->vertices())
      convex_to_fill_geometry->inputVertexPositions[v] = vec_to_GC_vec3(A * vec32vec(convex_to_fill_geometry->inputVertexPositions[v]));
    // Register the mesh with polyscope
    polyscope::SurfaceMesh *psHullFillMesh = polyscope::registerSurfaceMesh(
      "fillable hull",
      convex_to_fill_geometry->inputVertexPositions, convex_to_fill_mesh->getFaceVertexList(),
      polyscopePermutations(*convex_to_fill_mesh));
    psHullFillMesh->setTransparency(0.35);
      
    animate_convex_fill_deformation(mesh, geometry, convex_to_fill_mesh, convex_to_fill_geometry);
  }

  // IMPLOT attempts
  // bool open_flag = true;
  // ImGui::Begin("ImPlot Demo", &open_flag, ImGuiWindowFlags_MenuBar);
  // if (ImGui::BeginTabBar("ImPlotDemoTabs")) {
  //   if (ImGui::BeginTabItem("Plots")) {
  //     if (ImGui::TreeNodeEx("sub demo")) { 
  // ImGui::CreateContext();
  // ImGui::DestroyContext();       
  
  //       ImGui::TreePop();
  //     }
  //     ImGui::EndTabItem();
  //   }
  //   ImGui::EndTabBar();
  // }
  // ImGui::End();
  // DebugTextEncoding
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
  init_convex_shape_to_fill(all_polygons_current_item2);
  deformationSolver = new DeformationSolver(mesh, geometry, convex_to_fill_mesh, convex_to_fill_geometry);
  polyscope::state::userCallback = myCallback;

  // enrgy log stuff
  for(int i = 0; i < ENERGY_COUNT; i++)
        energies[i] = new float[MAX_FILL_ITERS];
  for(int i = 0; i < MAX_FILL_ITERS; i++)
    xs[i] = i;
  // polyscope::show();
  while (true) {
    if (animate){
      animate = false;
      gm_is_drawn = false;
      max_energy = 0.;
      energy_names[0] = "bending";
      energy_names[1] = "membrane";
      current_fill_iter = 0;
      current_ENERGY_COUNT = 3;
      Eigen::MatrixXd new_points = deformationSolver->solve_for_bending(1, true, &current_fill_iter, energies);
      // polyscope::SurfaceMesh *final_deformed_psMesh = polyscope::registerSurfaceMesh(
      //   "v2pipeline final mesh", new_points, mesh->getFaceVertexList(), polyscopePermutations(*mesh));

      // convex hull of deformed mesh obtained here
      Forward3DSolver* final_solver = new Forward3DSolver(mesh, new VertexPositionGeometry(*mesh, new_points), G, true);
      polyscope::SurfaceMesh *final_hull_psMesh = polyscope::registerSurfaceMesh(
        "v2pipeline final hull", final_solver->hullGeometry->inputVertexPositions, 
        final_solver->hullMesh->getFaceVertexList());
      final_hull_psMesh->setTransparency(0.4);

      // get the probability stuff
      final_solver->set_uniform_G();
      final_solver->updated = false;
      final_solver->initialize_pre_computes();
      BoundaryBuilder* tmp_bnd_builder = new BoundaryBuilder(final_solver);
      tmp_bnd_builder->build_boundary_normals();
      std::cout << "post-deform G:" << final_solver->get_G() << "\n";
      printf(" post deform G stuf:\n");
      update_visuals_with_G(final_solver, tmp_bnd_builder);
      tmp_bnd_builder->print_area_of_boundary_loops();
      printf("post deform dice energy: %f\n", tmp_bnd_builder->get_fair_dice_energy(fair_sides_count));
      
      
      if (pre_deform_G != Vector3({-10.,-10.,-10.})){// only if we coming here after dice energy search
        printf(" pre deform G stuff:\n");
        std::cout<< "pre-deform G:" << pre_deform_G << "\n";
        final_solver->set_G(pre_deform_G);
        final_solver->updated = false;
        final_solver->initialize_pre_computes();
        tmp_bnd_builder = new BoundaryBuilder(final_solver);
        tmp_bnd_builder->build_boundary_normals();
        vis_utils = VisualUtils(final_solver);
        update_visuals_with_G(final_solver, tmp_bnd_builder);
        // visualize_gauss_map();
        tmp_bnd_builder->print_area_of_boundary_loops();
        printf("pre deform G dice energy: %f\n", tmp_bnd_builder->get_fair_dice_energy(fair_sides_count));
        
      }
    }
    if (v2_dice_animate){
      max_energy = 0.01;
      current_fill_iter = 0;
      current_ENERGY_COUNT = 3;
      energy_names[0] = "bending";
      energy_names[1] = "membrane";
      // polyscope::removeAllStructures();
      version2_dice_pipeline(hull_opt_steps);
      v2_dice_animate = false;
    }
    polyscope::frameTick();
  }
  
  // convex_hull(forwardSolver->hullGeometry->inputVertexPositions);
  // build the solver
  
  // Give control to the polyscope gui
  

  return EXIT_SUCCESS;
}
