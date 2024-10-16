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
Vector3 G; // center of Mass

Forward3DSolver* forwardSolver;
  

VisualUtils vis_utils;

ManifoldSurfaceMesh* sphere_mesh;
VertexPositionGeometry* sphere_geometry;

// raster image stuff
FaceData<Vector3> face_colors;

// boundary stuff
BoundaryBuilder *boundary_builder;

bool compute_global_G_effect = true;
polyscope::PointCloud *test_pc;

int fair_sides_count = 6, // for optimization
    DE_step_count = 40;
bool do_sobolev_dice_grads = false,
     use_autodiff_for_dice_grad = true;
float sobolev_lambda = 2.,
      sobolev_lambda_decay = 0.95,
      dice_energy_step = 0.05,
      dice_search_decay = 0.98;
int sobolev_p = 2;
// optimization stuff



// example choice
std::vector<std::string> all_input_names = {std::string("triangular"), std::string("circus"), std::string("tet"), std::string("tet2"), std::string("cube")}; //{std::string("tet"), std::string("tet2"), std::string("cube"), std::string("tilted cube"), std::string("dodecahedron"), std::string("Conway spiral 4"), std::string("oloid")};
std::string input_name = "triangular";


void draw_stable_patches_on_gauss_map(bool on_height_surface = false, 
                                      BoundaryBuilder *bnd_builder = boundary_builder,
                                      Forward3DSolver *tmp_solver = forwardSolver,
                                      bool on_ambient_mesh = false){
  auto net_pair = build_and_draw_stable_patches_on_gauss_map(bnd_builder, vis_utils.center, vis_utils.gm_radi, vis_utils.arcs_seg_count, on_height_surface);
  if (on_ambient_mesh)
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
}


void init_visuals(){
  // Register the mesh with polyscope
  // psInputMesh = polyscope::registerSurfaceMesh(
  //   "init input mesh",
  //   geometry->inputVertexPositions, mesh->getFaceVertexList(),
  //   polyscopePermutations(*mesh));
  // psInputMesh->setTransparency(0.75);
  // psInputMesh->setEnabled(true);
  polyscope::SurfaceMesh *psHullMesh = polyscope::registerSurfaceMesh(
    "hull mesh",
    forwardSolver->hullGeometry->inputVertexPositions, forwardSolver->hullMesh->getFaceVertexList(),
    polyscopePermutations(*forwardSolver->hullMesh));
  forwardSolver->hullGeometry->requireFaceNormals();
  psHullMesh->addFaceVectorQuantity("face normals", forwardSolver->hullGeometry->faceNormals);
  psHullMesh->setEnabled(true);
  vis_utils.forwardSolver = forwardSolver;
  vis_utils.draw_G();
}

void update_solver(){ // only doing this for convex input
  forwardSolver = new Forward3DSolver(mesh, geometry, G, false);
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
  preprocess_mesh(mesh, geometry, triangulate, false, 1.);
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


std::vector<Eigen::Matrix3d> get_COM_grads_for_convex_uniform_shape(Eigen::MatrixX3d hull_positions){
  ManifoldSurfaceMesh *tmp_hull_mesh;
  VertexPositionGeometry *tmp_hull_geometry;
  std::tie(tmp_hull_mesh, tmp_hull_geometry) = get_mesh_for_convex_set(hull_positions);
  auto G_V_pair = find_center_of_mass(*tmp_hull_mesh, *tmp_hull_geometry);
  Vector3 tmp_G = G_V_pair.first;
  double volume = G_V_pair.second;

  Eigen::Matrix3d zmat = Eigen::Matrix3d::Zero();
  std::vector<Eigen::Matrix3d> dG_dv(hull_positions.rows());
  for (size_t i = 0; i < hull_positions.rows(); i++)
      dG_dv[i] = zmat;
  for (Face f: tmp_hull_mesh->faces()){
      // double face_area = tmp_solver->hullGeometry->faceArea(f);
      double face_area = polygonal_face_area(f, *tmp_hull_geometry);

      Vector3 face_normal = tmp_hull_geometry->faceNormal(f); // assuming outward normals
      size_t face_degree = f.degree();
      // assuming polygon faces here; 
      // TODO; check correctness for polygons
      Vector3 vert_sum = Vector3::zero();
      for (Vertex tmp_v: f.adjacentVertices())
          vert_sum += tmp_hull_geometry->inputVertexPositions[tmp_v];
      for (Halfedge he: f.adjacentHalfedges()){
          Vertex v = he.tailVertex();
          Vector3 p = tmp_hull_geometry->inputVertexPositions[v];
          Vector3 Gf_G = (vert_sum + p)/(double)(face_degree + 1) - tmp_G;
          DenseMatrix<double> tmp_mat = vec32vec(Gf_G) * 
                                        vec32vec(face_normal).transpose();
          assert(tmp_mat.cols() == 3);
          assert(tmp_mat.rows() == 3);
          dG_dv[v.getIndex()] += face_area * tmp_mat;
      }
  }
  for (size_t i = 0; i < hull_positions.rows(); i++)
      dG_dv[i] /= 3.*volume;
  return dG_dv;
}


void get_dice_energy_grads(Eigen::MatrixX3d hull_positions, Eigen::Vector3d G_vec,
                           Eigen::MatrixX3d &df_dv, Eigen::Vector3d &df_dG, double &dice_energy,
                           bool use_autodiff, bool frozen_G, std::string policy, FaceData<double> goal_probs, int fair_sides){
    
    if (!use_autodiff){
      Forward3DSolver fwd_solver(hull_positions, G_vec, true); // when getting grads, the input must be convex
      
      assert(hull_positions.rows() == fwd_solver.hullMesh->nVertices()); // check if input was convex

      if (!frozen_G){
        std::pair<Vector3, double> G_V_pair = find_center_of_mass(*fwd_solver.hullMesh, *fwd_solver.hullGeometry);
        fwd_solver.volume = G_V_pair.second;
        fwd_solver.set_G(G_V_pair.first); // since it deosn't make sense to gave non-uniform G with dependent G
      }
      
      fwd_solver.initialize_pre_computes();
      printf("   initialized for diffs\n");
      BoundaryBuilder *bnd_builder = new BoundaryBuilder(&fwd_solver);
      bnd_builder->build_boundary_normals();
      InverseSolver *inv_solver = new InverseSolver(bnd_builder);
      inv_solver->find_uni_mass_d_pf_dv(frozen_G);
      // distribution is set internally here

      // d_pf/dv
      inv_solver->set_fair_distribution_for_sink_faces(fair_sides); // top k faces are set
      VertexData<Vector3> approx_dice_energy_grads = inv_solver->find_uni_mass_total_vertex_grads(0.01);
      
      // d_pf/dG
      inv_solver->find_d_pf_d_Gs();
      Vector3 G_grad = inv_solver->find_total_g_grad();
      
      // convert to eigen
      dice_energy = bnd_builder->get_fair_dice_energy(fair_sides);

      // NOTE; forward solver reshuffles the vertices when taking hull, so need to map back
      for (Vertex v: fwd_solver.hullMesh->vertices()){
        df_dv.row(fwd_solver.org_hull_indices[v]) = vec32vec(approx_dice_energy_grads[v]);
      }
      df_dG = vec32vec(G_grad);
    }
    else { // autodiff
      // rebuild hull_positions since they were shuffled in solver's convex hull computation
      // Eigen::MatrixX3d hull_positions = vertex_data_to_matrix(fwd_solver.hullGeometry->inputVertexPositions); 
      auto dice_energy_lambda = [&] <typename Scalar> (const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &hull_poses_G_append_vec) -> Scalar {
        // decompose flat vector to positions and center of mass; G is the last 3 elements
        Eigen::Vector3<Scalar> G_eigen = hull_poses_G_append_vec.tail(3);
        size_t flat_n = hull_poses_G_append_vec.rows();
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 3> > hull_poses(hull_poses_G_append_vec.head(flat_n-3).data(), flat_n/3 - 1, 3);
        return BoundaryBuilder::dice_energy<Scalar>(hull_poses, G_eigen, policy, goal_probs, fair_sides, false);
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

      // populate df_dv by mapping to original input indices
      if (frozen_G){
        df_dv = dfdV;
      }
      else {
        std::vector<Eigen::Matrix3d> dG_dv = get_COM_grads_for_convex_uniform_shape(hull_positions);
        for (size_t i = 0; i < hull_positions.rows(); i++){
          df_dv.row(i) = dfdV.row(i) + (dG_dv[i].transpose() * df_dG).transpose();
        }
      }
    }
    
}


void dice_energy_opt(std::string policy, bool frozen_G, size_t step_count = 20){
  Forward3DSolver tmp_solver(mesh, geometry, G, true);
  tmp_solver.set_uniform_G();
  tmp_solver.initialize_pre_computes();
  
  G = tmp_solver.get_G();
  Eigen::Vector3d G_vec{G.x, G.y, G.z};
  Eigen::MatrixX3d hull_positions = vertex_data_to_matrix(tmp_solver.hullGeometry->inputVertexPositions);
  
  // BoundaryBuilder tmp_bnd_builder(&tmp_solver);
  
  double current_sobolev_lambda = sobolev_lambda;

  for (size_t iter = 0; iter < step_count; iter++){
    double dice_e;
    Eigen::Vector3d dfdG;
    Eigen::MatrixX3d dfdV(hull_positions.rows(), 3);
    
    printf("getting grads\n");
    FaceData<double> goal_probs;
    get_dice_energy_grads(hull_positions, G_vec, dfdV, dfdG, dice_e, 
                          use_autodiff_for_dice_grad, frozen_G,
                          policy, goal_probs, fair_sides_count);
    // update_visuals_with_G(&tmp_solver, &tmp_bnd_builder);
    std::cout << ANSI_FG_YELLOW << "i: "<< iter << "\tDE: " << dice_e << ANSI_RESET << "\n";
    polyscope::getPointCloud("Center of Mass")->updatePointPositions(std::vector<Vector3>{tmp_solver.get_G()});
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
    // polyscope::frameTick();
    // polyscope::screenshot(false);
    // polyscope::show();

    // printf("line search\n");
    double opt_step_size = hull_update_line_search(dfdV, hull_positions, G_vec, policy, goal_probs, fair_sides_count, 
                                                   dice_energy_step, dice_search_decay, frozen_G, 500);
    std::cout << ANSI_FG_RED << "  line search step size: " << opt_step_size << ANSI_RESET << "\n";
    hull_positions = hull_positions - opt_step_size * dfdV;
    tmp_solver = Forward3DSolver(hull_positions, G_vec, false); // could have been convaved
    if (!frozen_G){
      tmp_solver.set_uniform_G();
      G_vec = vec32vec(tmp_solver.get_G());
    }
    // IMPORTANT step; tmo solver's conv hull will shuffle the order of vertices
    hull_positions = vertex_data_to_matrix(tmp_solver.hullGeometry->inputVertexPositions);
    tmp_solver.initialize_pre_computes();
  }
  // Sorry; use boundary builder for visuals only

  BoundaryBuilder tmp_bnd_builder(&tmp_solver);
  tmp_bnd_builder.build_boundary_normals();
  update_visuals(&tmp_solver, &tmp_bnd_builder);
  tmp_bnd_builder.print_area_of_boundary_loops();



  // DEBUG for 11 sided double dice sum example
  FaceData<double> goal_probs = get_double_dice_probs_for_circus(&tmp_solver);
  FaceData<double> curr_probs_acum(*tmp_solver.hullMesh, 0.);
  for (Face f: tmp_solver.hullMesh->faces()){
    curr_probs_acum[f] = tmp_bnd_builder.face_region_area[tmp_solver.face_last_face[f]]/(4.*PI);
    if (goal_probs[f] > 0.){
      std::cout << "face " << f.getIndex() << "\tprob: " << tmp_bnd_builder.face_region_area[f]/(4.*PI) << "\tgoal:" << goal_probs[f] << " | " <<goal_probs[f]*36<<"/36 " << "\n";
    }
  }
  // polyscope::removeAllStructures();
  polyscope::registerSurfaceMesh("initial hull", geometry->inputVertexPositions, mesh->getFaceVertexList());
  polyscope::registerSurfaceMesh("current hull", tmp_solver.hullGeometry->inputVertexPositions, tmp_solver.hullMesh->getFaceVertexList());
  polyscope::registerPointCloud("Center of Mass", std::vector<Vector3>({tmp_solver.get_G()}))->setPointRadius(0.05);
  polyscope::getSurfaceMesh("current hull")->addFaceScalarQuantity(" current probs", tmp_bnd_builder.face_region_area/(4.*PI))->setEnabled(true);
  // writeSurfaceMesh(*tmp_solver.hullMesh, *tmp_solver.hullGeometry, "11_sided_double_dice_sum.obj");
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
  ImGui::SliderInt("ITERS", &DE_step_count, 1, 200);
  ImGui::SliderInt("fair sides", &fair_sides_count, 4, 20);
  ImGui::SliderFloat("DE step size", &dice_energy_step, 0, 0.5);
  ImGui::SliderFloat("DE step decay", &dice_search_decay, 0.1, 1.);
  
  if (ImGui::Button("dice energy opt")){
    dice_energy_opt("manual", false, DE_step_count);
  }
  if (ImGui::Button("show grads")){
    Forward3DSolver tmp_solver(mesh, geometry, G, true);
    tmp_solver.set_uniform_G();
    G = tmp_solver.get_G();
    tmp_solver.initialize_pre_computes();
    Eigen::Vector3d G_vec{G.x, G.y, G.z};
    Eigen::MatrixX3d hull_positions = vertex_data_to_matrix(tmp_solver.hullGeometry->inputVertexPositions);

    double dice_e;
    
    Eigen::Vector3d dfdG_ad, dfdG_proxy;
    Eigen::MatrixX3d dfdV_ad(hull_positions.rows(), 3), dfdV_proxy(hull_positions.rows(), 3);
    
    int fair_sides = fair_sides_count;

    printf(" geting unif G grads\n");
    bool frozen_G = false;
    std::string policy = "fair";
    FaceData<double> goal_probs;

    Eigen::Vector3d dfdG_ad_uniG, dfdG_proxy_uniG;
    Eigen::MatrixX3d dfdV_ad_uniG(hull_positions.rows(), 3), dfdV_proxy_uniG(hull_positions.rows(), 3);
    printf("getting AD grads\n");
    get_dice_energy_grads(hull_positions, G_vec, dfdV_ad_uniG, dfdG_ad_uniG, 
                          dice_e, true, frozen_G, policy, goal_probs, fair_sides);  
    printf("getting proxy grads\n");
    get_dice_energy_grads(hull_positions, G_vec, dfdV_proxy_uniG, dfdG_proxy_uniG, 
                          dice_e, false, frozen_G, policy, goal_probs, fair_sides);
    polyscope::getSurfaceMesh("hull mesh")->addVertexVectorQuantity("grads ad UNI", -1.*dfdV_ad_uniG)->setEnabled(true);
    polyscope::getSurfaceMesh("hull mesh")->addVertexVectorQuantity("grads proxy UNI", -1.*dfdV_proxy_uniG)->setEnabled(true);
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


// 