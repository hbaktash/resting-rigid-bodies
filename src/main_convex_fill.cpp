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
     frozen_G = false,
     structured_opt = true,
     one_time_CP_assignment = false,
     always_update_structure = true,
     with_hull_projection = false,
     use_initial_geometry = false,
     first_time = false;
bool do_remesh = false;
float scale_for_remesh = 2.;
polyscope::PointCloud *test_pc;

int fair_sides_count = 6; // for optimization


// optimization stuff
InverseSolver* inverseSolver;
float step_size = 0.01,
      step_size2 = 0.01,
      step_size3 = 0.5;
float stable_normal_update_thresh = -1;
int ARAP_max_iters = 10;

// deformation
DeformationSolver *deformationSolver;
bool animate = false,
     v2_dice_animate = false;
float scale_for_feasi = 2.;
float CP_lambda_exp = 1.,
      CP_mu = 1.15,
      barrier_lambda_exp = 1.,
      barrier_mu = 0.75;
int filling_max_iter = 100;
int hull_opt_steps = 100;
// example choice
std::vector<std::string> all_polyhedra_items = {std::string("tet"), std::string("tet2"), std::string("cube"), std::string("tilted cube"), std::string("sliced tet"), std::string("Conway spiral 4"), std::string("oloid"), std::string("fox"), std::string("small_bunny"), std::string("bunnylp"), std::string("kitten"), std::string("double-torus"), std::string("soccerball"), std::string("cowhead"), std::string("bunny"), std::string("gomboc")};
std::string all_polygons_current_item = "sliced tet",
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
    "init input mesh",
    geometry->inputVertexPositions, mesh->getFaceVertexList(),
    polyscopePermutations(*mesh));
  psInputMesh->setTransparency(0.75);
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
  convex_to_fill_geometry->inputVertexPositions *= scale_for_feasi;
  // convex_to_fill_geometry->inputVertexPositions += Vector3::constant(0.5);

  // Register the mesh with polyscope
  polyscope::SurfaceMesh *psHullFillMesh = polyscope::registerSurfaceMesh(
    "fillable hull",
    convex_to_fill_geometry->inputVertexPositions, convex_to_fill_mesh->getFaceVertexList(),
    polyscopePermutations(*convex_to_fill_mesh));
  psHullFillMesh->setTransparency(0.35);
}

void initialize_deformation_params(DeformationSolver *deformation_solver){
  deformation_solver->one_time_CP_assignment = one_time_CP_assignment;
  deformation_solver->CP_lambda = pow(10, CP_lambda_exp);
  deformation_solver->CP_mu = CP_mu;
  deformation_solver->barrier_init_lambda = pow(10, barrier_lambda_exp);
  deformation_solver->barrier_decay = barrier_mu;
  deformation_solver->filling_max_iter = filling_max_iter;  
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
         decay = 0.95,
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
          s *= decay;
  }
  printf("line search for dice ended at iter %d, s: %.10f, \n \t\t\t\t\t fnew: %f \n", j, s, tmp_fair_dice_energy);
  return found_smth_optimal ? s : 0.;
} 


void take_uni_mass_opt_vertices_step(bool frozen_G = false){
  if (!frozen_G){
    printf("finding uniform mass\n");
    forwardSolver->set_uniform_G();
  }
  printf("initalize precomputes\n");
  forwardSolver->initialize_pre_computes();
  printf("finding vertex derivatives\n");
  inverseSolver->find_uni_mass_d_pf_dv(frozen_G);
  printf(" vertex derivatives found!\n");
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
  double opt_step_size = hull_update_line_search(total_uni_mass_vertex_grads, forwardSolver, frozen_G);

  polyscope::registerSurfaceMesh("pre-step hull", forwardSolver->hullGeometry->inputVertexPositions,
                                                  forwardSolver->hullMesh->getFaceVertexList())->setEnabled(false);
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
      initialize_deformation_params(deformationSolver);

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
  if (!frozen_G){
    forwardSolver->set_uniform_G();
    G = forwardSolver->get_G();
  }
  update_solver_and_boundaries();
  update_visuals_with_G();
  boundary_builder->print_area_of_boundary_loops();
  printf(" fair dice energy: %f\n", boundary_builder->get_fair_dice_energy(fair_sides_count));
  // update_visuals_with_G();
}


void animate_convex_fill_deformation(ManifoldSurfaceMesh *_mesh, VertexPositionGeometry *_old_geometry,
                                     ManifoldSurfaceMesh *_convex_mesh, VertexPositionGeometry *_convex_geometry){
  deformationSolver = new DeformationSolver(_mesh, _old_geometry, _convex_mesh, _convex_geometry);   
  initialize_deformation_params(deformationSolver);
  animate = true;
}


void version2_dice_pipeline(size_t step_count = 1){
  // forwardSolver->set_uniform_G();
  ManifoldSurfaceMesh *hull_mesh = new ManifoldSurfaceMesh(forwardSolver->hullMesh->getFaceVertexList());
  VertexPositionGeometry *hull_geo = new VertexPositionGeometry(*hull_mesh, vertex_data_to_matrix(forwardSolver->hullGeometry->inputVertexPositions));
  Forward3DSolver *tmp_solver = new Forward3DSolver(hull_mesh, hull_geo, G);
  tmp_solver->set_uniform_G(); // frozen G matters here or not ???
  printf("initalize precomputes\n");
  tmp_solver->initialize_pre_computes();
  printf("finding vertex derivatives\n");
  
  polyscope::registerSurfaceMesh("pipe2 tmp sol", tmp_solver->inputGeometry->inputVertexPositions,
                                                  tmp_solver->inputMesh->getFaceVertexList());
  
  BoundaryBuilder *tmp_bnd_builder;
  InverseSolver *tmp_inv_solver;
  tmp_bnd_builder = new BoundaryBuilder(tmp_solver);
  tmp_inv_solver = new InverseSolver(tmp_bnd_builder);

  for (size_t iter = 0; iter < step_count; iter++){
    tmp_inv_solver->find_uni_mass_d_pf_dv(frozen_G);
    printf(" vertex derivatives found!\n");
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
    double opt_step_size = hull_update_line_search(vertex_grads, tmp_solver, frozen_G);
    
    auto [new_hull_mesh, new_hull_geo] = get_convex_hull_mesh(tmp_solver->hullGeometry->inputVertexPositions + opt_step_size * vertex_grads);
    tmp_solver = new Forward3DSolver(new_hull_mesh, new_hull_geo, tmp_solver->get_G(), false); 
    // TODO: is this ok?
  
    // visuals
    polyscope::registerSurfaceMesh("pipe2 tmp sol", tmp_solver->inputGeometry->inputVertexPositions,
                                                  tmp_solver->inputMesh->getFaceVertexList());
    polyscope::frameTick();


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


  polyscope::warning(" dice enrgy search done!", " deforming initial shape");
  polyscope::getSurfaceMesh("pipe2 tmp sol")->setEnabled(false);
  polyscope::getSurfaceMesh("init input mesh")->setEnabled(false);
  polyscope::getSurfaceMesh("init hull mesh")->setEnabled(false);


  // Register the mesh with polyscope
  polyscope::SurfaceMesh *psHullFillMesh = polyscope::registerSurfaceMesh(
    "fillable hull", tmp_solver->inputGeometry->inputVertexPositions, tmp_solver->inputMesh->getFaceVertexList());
  psHullFillMesh->setTransparency(0.35);

  deformationSolver = new DeformationSolver(mesh, geometry, tmp_solver->inputMesh, tmp_solver->inputGeometry);   
  initialize_deformation_params(deformationSolver);
  DenseMatrix<double> new_points = deformationSolver->solve_for_bending(1);

  polyscope::SurfaceMesh *final_deformed_psMesh = polyscope::registerSurfaceMesh(
    "v2pipeline final mesh", new_points, mesh->getFaceVertexList(), polyscopePermutations(*mesh));

  std::cout << "pre-deform G:" << tmp_solver->get_G() << "\n";
  // convex hull of deformed mesh obtained here
  Forward3DSolver* final_solver = new Forward3DSolver(mesh, new VertexPositionGeometry(*mesh, new_points), G, true);
  polyscope::SurfaceMesh *final_hull_psMesh = polyscope::registerSurfaceMesh(
    "v2pipeline final hull", final_solver->hullGeometry->inputVertexPositions, 
    final_solver->hullMesh->getFaceVertexList());
  final_hull_psMesh->setTransparency(0.4);

  // get the probability stuff
  final_solver->set_uniform_G();
  std::cout << "post-deform G:" << final_solver->get_G() << "\n";
  final_solver->updated = false;
  final_solver->initialize_pre_computes();
  tmp_bnd_builder = new BoundaryBuilder(final_solver);
  tmp_bnd_builder->build_boundary_normals();
  update_visuals_with_G(final_solver, tmp_bnd_builder);
  printf(" post deform probabilities:\n");
  tmp_bnd_builder->print_area_of_boundary_loops();
  printf("post deform dice energy: %f\n", tmp_bnd_builder->get_fair_dice_energy(fair_sides_count));
  
  polyscope::warning("testing G effect\n");
  final_solver->set_G(tmp_solver->get_G());
  final_solver->updated = false;
  final_solver->initialize_pre_computes();
  tmp_bnd_builder = new BoundaryBuilder(final_solver);
  tmp_bnd_builder->build_boundary_normals();
  printf(" good G probabilities:\n");
  tmp_bnd_builder->print_area_of_boundary_loops();
  printf("good G dice energy: %f\n", tmp_bnd_builder->get_fair_dice_energy(fair_sides_count));
  
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
          }
          if (is_selected)
              ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
      }
      ImGui::EndCombo();
  }

  if (ImGui::SliderFloat("scale for feasibility", &scale_for_feasi, 1., 10.));
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
    
    animate_convex_fill_deformation(mesh, geometry, convex_to_fill_mesh, convex_to_fill_geometry);
    // TODO: check what should be done here to avoid the ugly bool trick
    // deformationSolver->solve_for_bending(1);
  }
  if (ImGui::Checkbox("one time CP assignment", &one_time_CP_assignment)) deformationSolver->one_time_CP_assignment = one_time_CP_assignment;
  if (ImGui::SliderFloat("CP lambda log10/initial value", &CP_lambda_exp, 0., 5.)) deformationSolver->CP_lambda = pow(10, CP_lambda_exp);
  if (ImGui::SliderFloat("CP growth mu ", &CP_mu, 1., 1.5)) deformationSolver->CP_mu = CP_mu;
  if (ImGui::SliderFloat("Barrier lambda log10/initial value", &barrier_lambda_exp, 0., 5.)) deformationSolver->barrier_init_lambda = pow(10, barrier_lambda_exp);
  if (ImGui::SliderFloat("Barrier decay mu ", &barrier_mu, 0., 1.)) deformationSolver->barrier_decay = barrier_mu;
  
  if (ImGui::SliderInt("filling iters", &filling_max_iter, 0, 200)) deformationSolver->filling_max_iter = filling_max_iter;

  if (ImGui::Button("uniform mass G")){
    G = find_center_of_mass(*forwardSolver->inputMesh, *forwardSolver->inputGeometry).first;
    update_solver_and_boundaries();
    update_visuals_with_G();
  }
  if (ImGui::Checkbox("draw artificial R3 boundaries", &test_guess)) draw_stable_patches_on_gauss_map();
  
  if (ImGui::Button("take fair step (joint gradient)")) {
    take_uni_mass_opt_vertices_step(frozen_G);
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

  if (ImGui::Button("optimize hull (w.r.t. dice energy)")) {
    v2_dice_animate = true;
    if (polyscope::hasSurfaceMesh("v2pipeline final mesh")) polyscope::removeSurfaceMesh("v2pipeline final mesh");
    if (polyscope::hasSurfaceMesh("v2pipeline final hull")) polyscope::removeSurfaceMesh("v2pipeline final hull");
  }
  if (ImGui::SliderInt("hull optimize step count", &hull_opt_steps, 1, 200));
  
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

  // polyscope::show();
  while (true) {
    if (animate){
      deformationSolver->solve_for_bending(1);
      animate = false;
    }
    if (v2_dice_animate){
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
