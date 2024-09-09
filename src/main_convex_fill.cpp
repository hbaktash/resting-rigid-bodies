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
#include "markov_model.h"
#include "inv_design.h"
// #include "igl/arap.h"
#include "deformation.h"
#include "implot.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


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
     first_time = false,
     curvature_weighted_CP = false;
bool do_remesh = false;
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
int ARAP_max_iters = 10;

// deformation
DeformationSolver *deformationSolver;
bool animate = false,
     animate_G_deform = false,
     v2_dice_animate = false,
     enforce_snapping = false;
float dice_search_decay = 0.95;

float bending_lambda_exps[2] = {1., 1.},
      membrane_lambda_exps[2] = {3., 3.},
      CP_lambda_exps[2] = {1., 7.},
      barrier_lambda_exps[2] = {-4., -8.},
      G_lambda_exps[2] = {1,5},
      reg_lambda_exp = -3.,
      internal_p = 0.91,
      refinement_CP_threshold = 0.001,
      active_set_threshold = 0.08,
      split_robustness_threshold = 0.2;
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
    auto net_pair = build_and_draw_stable_patches_on_gauss_map(bnd_builder, vis_utils.center, vis_utils.gm_radi, vis_utils.arcs_seg_count, on_height_surface);
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
  deformation_solver->active_set_threshold = active_set_threshold;
  deformation_solver->refinement_CP_threshold = refinement_CP_threshold;
  deformation_solver->split_robustness_threshold = split_robustness_threshold;
  deformation_solver->enforce_snapping = enforce_snapping;
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
    std::vector<Face> stable_faces;
    for (Face f: fwd_solver->hullMesh->faces()){
      if (fwd_solver->face_is_stable(f))
        stable_faces.push_back(f);
    }
    face_colors = generate_random_colors(fwd_solver->hullMesh, stable_faces);
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


void animate_convex_fill_deformation(ManifoldSurfaceMesh *_mesh, VertexPositionGeometry *_old_geometry,
                                     ManifoldSurfaceMesh *_convex_mesh, VertexPositionGeometry *_convex_geometry){
  if (polyscope::hasSurfaceMesh("fillable hull")) polyscope::getSurfaceMesh("fillable hull")->setEnabled(true);
  deformationSolver = new DeformationSolver(_mesh, _old_geometry, _convex_mesh, _convex_geometry);   
  initialize_deformation_params(deformationSolver);
  animate = true;
}

void animate_G_diff_deformation(Vector3 ideal_G){
  if (polyscope::hasSurfaceMesh("fillable hull")) polyscope::getSurfaceMesh("fillable hull")->setEnabled(true);
  deformationSolver->goal_G = ideal_G;   
  deformationSolver->init_G_lambda = pow(10, G_lambda_exps[0]);
  deformationSolver->final_G_lambda = pow(10, G_lambda_exps[1]);
  initialize_deformation_params(deformationSolver);
  animate_G_deform = true;
}

void dice_energy_opt_stan(size_t step_count = 1){
  // forwardSolver->set_uniform_G();
  // Eigen::MatrixX3d point_cloud = vertex_data_to_matrix(points);
  Forward3DSolver tmp_solver(mesh, geometry, G, true);
  tmp_solver.set_uniform_G();
  G = tmp_solver.get_G();
  Eigen::Vector3d G_vec{G.x, G.y, G.z};
  Eigen::MatrixX3d point_cloud = vertex_data_to_matrix(tmp_solver.hullGeometry->inputVertexPositions);
  std::cout << "initial G: " << G_vec << "\n";
  tmp_solver.initialize_pre_computes();
  // if (polyscope::hasSurfaceMesh("fillable hull")) polyscope::getSurfaceMesh("fillable hull")->setEnabled(false);
  // if (polyscope::hasSurfaceMesh("temp sol")) polyscope::getSurfaceMesh("temp sol")->setEnabled(false);
  // if (polyscope::hasSurfaceMesh("init hull mesh")) polyscope::getSurfaceMesh("init hull mesh")->setEnabled(false);
  // if (polyscope::hasSurfaceMesh("init input mesh")) polyscope::getSurfaceMesh("init input mesh")->setEnabled(false);
  // polyscope::screenshot(false);
  BoundaryBuilder tmp_bnd_builder(&tmp_solver);
  
  max_energy = 0.;
  current_fill_iter = -1;
  current_ENERGY_COUNT = 3;
  energy_names[0] = "Sblv grad norm";
  energy_names[1] = "raw grad norm";
  energy_names[2] = "Sblv/raw grad norm";
  double current_sobolev_lambda = sobolev_lambda;

  for (size_t iter = 0; iter < step_count; iter++){
    tmp_bnd_builder.build_boundary_normals();
    auto dice_energy_lambda = [&] <typename Scalar> (const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &hull_poses_G_append_vec) -> Scalar {
      // decompose flat vector to positions and center of mass; G is the last 3 elements
      Eigen::Vector3<Scalar> G_eigen = hull_poses_G_append_vec.tail(3);
      size_t flat_n = hull_poses_G_append_vec.rows();
      Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 3> > hull_poses(hull_poses_G_append_vec.head(flat_n-3).data(), flat_n/3 - 1, 3);
      return BoundaryBuilder::dice_energy<Scalar>(hull_poses, G_eigen, tmp_solver, fair_sides_count);
    };
    Eigen::VectorXd hull_poses_vec = point_cloud.reshaped();
    Eigen::VectorXd hull_poses_and_G_vec(hull_poses_vec.size() + 3);
    hull_poses_and_G_vec << hull_poses_vec, G_vec;
    
    Eigen::VectorXd dfdU_vec;
    double dice_e;
    stan::math::gradient(dice_energy_lambda, hull_poses_and_G_vec, dice_e, dfdU_vec);
    Eigen::Vector3d G_grad = dfdU_vec.tail(3);
    size_t flat_n = dfdU_vec.rows();
    Eigen::Map<Eigen::MatrixXd> dfdV(dfdU_vec.head(flat_n-3).data(), flat_n/3 - 1, 3);

    // update_visuals_with_G(&tmp_solver, &tmp_bnd_builder);
    std::cout << ANSI_FG_YELLOW << " DE at iter "<< iter << " f: " << dice_e << ANSI_RESET << "\n";
    auto curr_hull_psmesh = polyscope::registerSurfaceMesh("current hull", 
                                                  tmp_solver.hullGeometry->inputVertexPositions,
                                                  tmp_solver.hullMesh->getFaceVertexList());
    curr_hull_psmesh->setEnabled(true);
    curr_hull_psmesh->setEdgeWidth(2.);
    curr_hull_psmesh->setTransparency(0.7);

    auto old_tmp_geo = tmp_solver.inputGeometry->inputVertexPositions;
    auto old_tmp_mesh = tmp_solver.inputMesh->getFaceVertexList();
    
    // show grads
    polyscope::getSurfaceMesh("current hull")->addVertexVectorQuantity("ad total grads", dfdV*(-1.))->setEnabled(true);
    if (do_sobolev_dice_grads){
      
      current_fill_iter = iter;
      current_sobolev_lambda *= sobolev_lambda_decay;
      Eigen::MatrixXd diffused_dfdV = sobolev_diffuse_gradients(dfdV, *tmp_solver.hullMesh, sobolev_lambda, sobolev_p);
      
      double raw_norm = tinyAD_flatten(dfdV).norm(); // use eigen flaten
      double diffused_norm = tinyAD_flatten(diffused_dfdV).norm();
      energies[0][iter] = tmp_bnd_builder.get_fair_dice_energy(fair_sides_count);
      energies[1][iter] = diffused_norm / raw_norm;
      
      // std::cout << "energies: " << energies[0][iter] << " " << energies[1][iter] << " " << energies[2][iter] << "\n";
      curr_hull_psmesh->addVertexVectorQuantity(" sobolev diffused grads", diffused_dfdV*(-1.))->setEnabled(true);
      dfdV = diffused_dfdV;
    }

    curr_hull_psmesh->setSurfaceColor({0.1,0.9,0.1});
    curr_hull_psmesh->setEdgeWidth(2.);
    dfdV = dfdV * -1.;
    double opt_step_size = hull_update_line_search(vertex_matrix_to_data(dfdV, *tmp_solver.hullMesh), tmp_solver, fair_sides_count, 
                                                   step_size3, dice_search_decay, frozen_G);
    std::cout << ANSI_FG_RED << " sobolev lambda: " << current_sobolev_lambda 
              << " \n opt step size: " << opt_step_size << ANSI_RESET << "\n";
    auto [new_hull_mesh, new_hull_geo] = get_convex_hull_mesh(point_cloud + opt_step_size * dfdV);
    tmp_solver = Forward3DSolver(new_hull_mesh, new_hull_geo, tmp_solver.get_G(), false); 

    //DEBUG
    polyscope::frameTick();
    // polyscope::screenshot(false);
    // polyscope::show();
    
    if (!frozen_G){
      tmp_solver.set_uniform_G();
      Vector3 GG = tmp_solver.get_G();
      G_vec = {GG.x, GG.y, GG.z};
    }
    tmp_solver.updated = false;
    tmp_solver.initialize_pre_computes();
    tmp_bnd_builder = BoundaryBuilder(&tmp_solver);

    // update vertorized positions
    point_cloud = vertex_data_to_matrix(tmp_solver.inputGeometry->inputVertexPositions);
  }
  tmp_bnd_builder.build_boundary_normals();
  update_visuals_with_G(&tmp_solver, &tmp_bnd_builder);
  tmp_bnd_builder.print_area_of_boundary_loops();
  draw_boundary_patches = true;
  test_guess = true;
  draw_stable_patches_on_gauss_map(false, &tmp_bnd_builder);
  // //DEBUG
  // if (polyscope::hasSurfaceMesh("init hull mesh")) polyscope::getSurfaceMesh("init hull mesh")->setEnabled(true);

  // polyscope::removeAllStructures();
  // polyscope::SurfaceMesh *psHullFillMesh = polyscope::registerSurfaceMesh(
  //   "fillable hull", tmp_solver.inputGeometry->inputVertexPositions, tmp_solver.inputMesh->getFaceVertexList());
  // psHullFillMesh->setTransparency(0.45);
  // psHullFillMesh->setEdgeWidth(2.);
  // psHullFillMesh->setEnabled(true);
  // gm_is_drawn = false;
  convex_to_fill_mesh = tmp_solver.hullMesh;
  convex_to_fill_geometry = tmp_solver.hullGeometry;
  // load_hull_from_file = false;

  // pre_deform_G = tmp_solver.get_G();
  // std::cout << "pre-deform G:" << tmp_solver.get_G() << "\n";
}


void dice_energy_opt_proxy(size_t step_count = 1){
  // forwardSolver->set_uniform_G();
  Forward3DSolver tmp_solver(mesh, geometry, G, true);
  tmp_solver.set_uniform_G();
  G = tmp_solver.get_G();
  Eigen::Vector3d G_vec{G.x, G.y, G.z};
  Eigen::MatrixX3d point_cloud = vertex_data_to_matrix(tmp_solver.hullGeometry->inputVertexPositions);
  std::cout << "initial G: " << G_vec << "\n";
  tmp_solver.initialize_pre_computes();
  BoundaryBuilder tmp_bnd_builder(&tmp_solver);
  InverseSolver tmp_inv_solver(&tmp_bnd_builder);
  
  // if (polyscope::hasSurfaceMesh("fillable hull")) polyscope::getSurfaceMesh("fillable hull")->setEnabled(false);
  // if (polyscope::hasSurfaceMesh("temp sol")) polyscope::getSurfaceMesh("temp sol")->setEnabled(false);
  // if (polyscope::hasSurfaceMesh("init hull mesh")) polyscope::getSurfaceMesh("init hull mesh")->setEnabled(false);
  // if (polyscope::hasSurfaceMesh("init input mesh")) polyscope::getSurfaceMesh("init input mesh")->setEnabled(false);
  // auto pip2psmesh = polyscope::registerSurfaceMesh("pipe2 tmp sol", tmp_solver->inputGeometry->inputVertexPositions,
  //                                                 tmp_solver->inputMesh->getFaceVertexList());
  // pip2psmesh->setEdgeWidth(2.);
  // polyscope::screenshot(false);
  // BoundaryBuilder tmp_bnd_builder(&tmp_solver);

  max_energy = 0.;
  // cur = 2
  current_fill_iter = -1;
  current_ENERGY_COUNT = 2;
  energy_names[0] = "Dice energy";
  energy_names[1] = "Sblv/raw grad norm";
  double current_sobolev_lambda = sobolev_lambda;
  for (size_t iter = 0; iter < step_count; iter++){

    printf("finding vertex derivatives\n");
    tmp_bnd_builder.build_boundary_normals(); // face-last-face is called 
    
    update_visuals_with_G(&tmp_solver, &tmp_bnd_builder); // TODO
    std::cout << ANSI_FG_YELLOW << " fair dice energy iter "<< iter<< "  f: " << tmp_bnd_builder.get_fair_dice_energy(fair_sides_count) << ANSI_RESET << "\n";

    auto pip2psmesh = polyscope::registerSurfaceMesh("pipe2 tmp sol", tmp_solver.hullGeometry->inputVertexPositions,
                                                  tmp_solver.hullMesh->getFaceVertexList());

    if (first_time || always_update_structure){
      if (first_time) stable_normal_update_thresh = -1.;
      else stable_normal_update_thresh = 0.1;
      first_time = false;
    }

    // checking approx grad vs real grads
    tmp_inv_solver.find_uni_mass_d_pf_dv(false, frozen_G);
    VertexData<Vector3> approx_dice_energy_grads = tmp_inv_solver.find_uni_mass_total_vertex_grads(fair_sides_count,
                                                                                                   stable_normal_update_thresh);
    
    // tmp_inv_solver->find_uni_mass_d_pf_dv(use_autodiff_for_dice_grad, frozen_G);
    // VertexData<Vector3> dice_energy_grads = tmp_inv_solver->find_uni_mass_total_vertex_grads(fair_sides_count,
    //                                                                                     structured_opt, stable_normal_update_thresh);
    auto old_tmp_geo = tmp_solver.inputGeometry->inputVertexPositions;
    auto old_tmp_mesh = tmp_solver.inputMesh->getFaceVertexList();
    
    pip2psmesh->setEnabled(true);
    pip2psmesh->setTransparency(0.7);
    // DEBUG
    pip2psmesh->addVertexVectorQuantity("approx total grads", approx_dice_energy_grads)->setEnabled(true);
    // polyscope::show();

    if (do_sobolev_dice_grads){
      current_fill_iter = iter;
      double raw_norm = tinyAD_flatten(vertex_data_to_matrix(approx_dice_energy_grads)).norm();
      current_sobolev_lambda *= sobolev_lambda_decay;
      approx_dice_energy_grads = sobolev_diffuse_gradients(approx_dice_energy_grads, *tmp_solver.hullMesh, current_sobolev_lambda, sobolev_p);
      // dice_energy_grads = tmp_inv_solver->sobolev_diffuse_gradients(dice_energy_grads, current_sobolev_lambda, sobolev_p);
      
      double sobolev_norm = tinyAD_flatten(vertex_data_to_matrix(approx_dice_energy_grads)).norm();
      
      energies[0][iter] = tmp_bnd_builder.get_fair_dice_energy(fair_sides_count);
      energies[1][iter] = sobolev_norm / raw_norm;
      // std::cout << "energies: " << energies[0][iter] << " " << energies[1][iter] << " " << energies[2][iter] << "\n";
      pip2psmesh->addVertexVectorQuantity(" sobolev diffused dice grads", approx_dice_energy_grads)->setEnabled(true);
    }

    pip2psmesh->setSurfaceColor({0.1,0.9,0.1});
    pip2psmesh->setEdgeWidth(2.);
    // polyscope::screenshot(false);
    polyscope::frameTick();
    // polyscope::show();

    double opt_step_size = hull_update_line_search(approx_dice_energy_grads, tmp_solver, fair_sides_count, step_size3, dice_search_decay, frozen_G);
    std::cout << ANSI_FG_RED << " sobolev lambda: " << current_sobolev_lambda 
              << "opt step size: " << opt_step_size << ANSI_RESET << "\n";
    auto [new_hull_mesh, new_hull_geo] = get_convex_hull_mesh(tmp_solver.hullGeometry->inputVertexPositions + opt_step_size * approx_dice_energy_grads);
    tmp_solver = Forward3DSolver(new_hull_mesh, new_hull_geo, tmp_solver.get_G(), false); 
    // TODO: is this ok?
  
    // visuals
    if (!frozen_G){
      tmp_solver.set_uniform_G();
    }
    tmp_solver.updated = false;
    tmp_solver.initialize_pre_computes();
    tmp_bnd_builder = BoundaryBuilder(&tmp_solver);
    tmp_inv_solver  = InverseSolver(&tmp_bnd_builder);
  }
  tmp_bnd_builder.build_boundary_normals();
  update_visuals_with_G(&tmp_solver, &tmp_bnd_builder);
  tmp_bnd_builder.print_area_of_boundary_loops();
  draw_boundary_patches = true;
  test_guess = true;
  draw_stable_patches_on_gauss_map(false, &tmp_bnd_builder);
  convex_to_fill_mesh = tmp_solver.hullMesh;
  convex_to_fill_geometry = tmp_solver.hullGeometry;// gm_is_drawn = false;
  
  // //DEBUG
  // if (polyscope::hasSurfaceMesh("init hull mesh")) polyscope::getSurfaceMesh("init hull mesh")->setEnabled(true);

  // polyscope::removeAllStructures();
  // polyscope::SurfaceMesh *psHullFillMesh = polyscope::registerSurfaceMesh(
  //   "fillable hull", tmp_solver.inputGeometry->inputVertexPositions, tmp_solver.inputMesh->getFaceVertexList());
  // psHullFillMesh->setTransparency(0.45);
  // psHullFillMesh->setEdgeWidth(2.);
  // psHullFillMesh->setEnabled(true);
  // load_hull_from_file = false;
  
  // pre_deform_G = tmp_solver->get_G();
  // std::cout << "pre-deform G:" << tmp_solver->get_G() << "\n";
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
  if (ImGui::Button("deform for G differece")){
    polyscope::removeAllStructures();

    polyscope::SurfaceMesh *psHullFillMesh = polyscope::registerSurfaceMesh(
    "fillable hull",
    convex_to_fill_geometry->inputVertexPositions, convex_to_fill_mesh->getFaceVertexList(),
    polyscopePermutations(*convex_to_fill_mesh));
    psHullFillMesh->setTransparency(0.35);
    animate_G_diff_deformation(pre_deform_G);
  }
  if (ImGui::InputFloat2("init/final Gdiff log ", G_lambda_exps));
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
  if (ImGui::Checkbox("enforce snapping at threshold", &enforce_snapping));
  if (ImGui::SliderFloat("split rubostness threshold ", &split_robustness_threshold, 0., 1.));
  if (ImGui::InputFloat2("init/final bending log ", bending_lambda_exps)
    || ImGui::InputFloat2("init/final membrane log ", membrane_lambda_exps)
    || ImGui::InputFloat2("init/final CP log ", CP_lambda_exps)
    || ImGui::Checkbox("use QP solver ", &use_QP_solver));
  if(!use_QP_solver) 
    if (ImGui::InputFloat2("init/final barrier log ", barrier_lambda_exps));
  else
    if (ImGui::InputFloat("active set threshold ", &active_set_threshold));
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
  if (ImGui::Checkbox("frozen G", &frozen_G));

  // // V3 stuff 
  // if (ImGui::Button("take fair step (joint gradient)")) {
  //   take_uni_mass_opt_vertices_step(frozen_G);
  // }
  // if (ImGui::Checkbox("deform after", &deform_after));
  // if (ImGui::Checkbox("compute global G effect", &compute_global_G_effect)) inverseSolver->compute_global_G_effect = compute_global_G_effect;
  // if (ImGui::Checkbox("structured opt", &structured_opt));
  // if (ImGui::Checkbox("update structured at every step", &always_update_structure));
  // if (ImGui::Checkbox("decrease convex hull points", &with_hull_projection));

  if (ImGui::SliderFloat("DE line search decay", &dice_search_decay, 0., 1.));
  if (ImGui::SliderInt("DE iterations", &hull_opt_steps, 1, 200));
  if (ImGui::SliderFloat("stable normal update thresh", &stable_normal_update_thresh, 0., 4.0));
  if (ImGui::Checkbox("Sobolev pre-condition", &do_sobolev_dice_grads));
  if (do_sobolev_dice_grads){
    if (ImGui::SliderFloat("Sobolev lambda", &sobolev_lambda, 0., 5.));
    if (ImGui::SliderFloat("Sobolev lambda decay", &sobolev_lambda_decay, 0., 1));
    if (ImGui::SliderInt("Sobolev power", &sobolev_p, 1, 5));
  }
  if (ImGui::Checkbox("autodiff DE grads", &use_autodiff_for_dice_grad));
  if (ImGui::Button("optimize hull DE")) {
    v2_dice_animate = true;
    if (polyscope::hasSurfaceMesh("v2pipeline final mesh")) polyscope::removeSurfaceMesh("v2pipeline final mesh");
    if (polyscope::hasSurfaceMesh("v2pipeline final hull")) polyscope::removeSurfaceMesh("v2pipeline final hull");
  }
  if (ImGui::Button("Save hull and G to file")){
    writeSurfaceMesh(*convex_to_fill_mesh, *convex_to_fill_geometry, "../meshes/hulls/hull_" + all_polygons_current_item + "_fs" + std::to_string(fair_sides_count), "obj");
    std::string filename = "G_" + all_polygons_current_item + "_fs" + std::to_string(fair_sides_count) + ".txt";
    std::ofstream file(filename);
    file << G.x << " " << G.y << " " << G.z << "\n";
    file.close();
  }

  // Test stuff
  // if ( ImGui::InputFloat3("row0", row0)
  //   || ImGui::InputFloat3("row1", row1)
  //   || ImGui::InputFloat3("row2", row2));
  // if (ImGui::Button("simple test")){
  //   DenseMatrix<double> A = DenseMatrix<double>::Zero(3,3);
  //   A(0,0) = row0[0]; A(0,1) = row0[1]; A(0,2) = row0[2];
  //   A(1,0) = row1[0]; A(1,1) = row1[1]; A(1,2) = row1[2];
  //   A(2,0) = row2[0]; A(2,1) = row2[1]; A(2,2) = row2[2];
  //   std::cout << "A matrix\n \t" << A << "\n";
  //   // auto tmp_def = new DeformationSolver(forwardSolver->inputMesh, forwardSolver->inputGeometry, convex_to_fill_mesh, convex_to_fill_geometry);   
  //   // initialize_deformation_params(tmp_def);
  //   // tmp_def->print_energies_after_transform(A);
  //   polyscope::removeAllStructures();
  //   convex_to_fill_mesh = forwardSolver->hullMesh;
  //   convex_to_fill_geometry = forwardSolver->hullGeometry;
  //   preprocess_mesh(convex_to_fill_mesh, convex_to_fill_geometry, true);

  //   // transform
  //   for (Vertex v: convex_to_fill_mesh->vertices())
  //     convex_to_fill_geometry->inputVertexPositions[v] = vec_to_GC_vec3(A * vec32vec(convex_to_fill_geometry->inputVertexPositions[v]));
  //   // Register the mesh with polyscope
  //   polyscope::SurfaceMesh *psHullFillMesh = polyscope::registerSurfaceMesh(
  //     "fillable hull",
  //     convex_to_fill_geometry->inputVertexPositions, convex_to_fill_mesh->getFaceVertexList(),
  //     polyscopePermutations(*convex_to_fill_mesh));
  //   psHullFillMesh->setTransparency(0.35);
      
  //   animate_convex_fill_deformation(mesh, geometry, convex_to_fill_mesh, convex_to_fill_geometry);
  // }
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
      energy_names[2] = "CP";
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
      post_deform_G = final_solver->get_G();
      std::cout << "post-deform G:" << post_deform_G << "\n";
      printf(" post deform G stuf:\n");
      update_visuals_with_G(final_solver, tmp_bnd_builder);
      tmp_bnd_builder->print_area_of_boundary_loops();
      printf("post deform dice energy: %f\n", tmp_bnd_builder->get_fair_dice_energy(fair_sides_count));
      
      bool test_predeformation = true;
      if (pre_deform_G == Vector3({-10.,-10.,-10.})){
        test_predeformation = false;
        pre_deform_G = find_center_of_mass(*deformationSolver->convex_mesh, *deformationSolver->convex_geometry).first;
      }
      if (test_predeformation){
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
    if (animate_G_deform){
      animate_G_deform = false;
      gm_is_drawn = false;
      max_energy = 0.1;
      energy_names[0] = "bending";
      energy_names[1] = "membrane";
      energy_names[2] = "G diff";
      current_fill_iter = 0;
      current_ENERGY_COUNT = 3;
      Eigen::MatrixXd new_points = deformationSolver->solve_for_G(1, true, &current_fill_iter, energies);

      Forward3DSolver* final_solver = new Forward3DSolver(mesh, new VertexPositionGeometry(*mesh, new_points), G, true);
      // get the probability stuff
      final_solver->set_uniform_G();
      final_solver->updated = false;
      final_solver->initialize_pre_computes();
      BoundaryBuilder* tmp_bnd_builder = new BoundaryBuilder(final_solver);
      tmp_bnd_builder->build_boundary_normals();
      printf(" post G fix stuf:\n");
      update_visuals_with_G(final_solver, tmp_bnd_builder);
      tmp_bnd_builder->print_area_of_boundary_loops();
      printf("post G fix dice energy: %f\n", tmp_bnd_builder->get_fair_dice_energy(fair_sides_count));
    }
    if (v2_dice_animate){
      max_energy = 0.01;
      current_fill_iter = 0;
      current_ENERGY_COUNT = 3;
      energy_names[0] = "bending";
      energy_names[1] = "membrane";
      // polyscope::removeAllStructures();
      if (use_autodiff_for_dice_grad){
        forwardSolver->set_uniform_G();
        dice_energy_opt_stan(hull_opt_steps);
      }
      else
        dice_energy_opt_proxy(hull_opt_steps);
      v2_dice_animate = false;
    }
    polyscope::frameTick();
  }

  return EXIT_SUCCESS;
}
