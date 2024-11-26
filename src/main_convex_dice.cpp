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
std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr;
ManifoldSurfaceMesh *mesh, *optimized_mesh;
VertexPositionGeometry *geometry, *optimized_geometry;
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
bool do_sobolev_dice_grads = true,
     use_autodiff_for_dice_grad = true,
     visualize_steps = true,
     adaptive_reg = false,
     frozen_G = true;
float sobolev_lambda = 2.,
      sobolev_lambda_decay = 0.8,
      dice_energy_step = 0.01,
      dice_search_decay = 0.98,
      bary_reg = 0.1,
      coplanar_reg = 0.0,
      cluster_distance_reg = 0.0;
int sobolev_p = 2;
// optimization stuff



// example choice
std::vector<std::string> all_input_names = {std::string("6 prism"), std::string("hendecahedron"), std::string("triangular"), std::string("circus"), std::string("icosahedron"), std::string("dodecahedron"), std::string("cuub"), std::string("octahedron")}; // {std::string("tet"), std::string("tet2"), std::string("cube"), std::string("tilted cube"), std::string("dodecahedron"), std::string("Conway spiral 4"), std::string("oloid")};
std::string input_name = "7 prism";
std::string policy_general = "manualCluster"; // "fair", "manualCluster ", "manual"
// std::string policy_shape = "dodecahedron binomial"; // "dodecahedron binomial", "octahedron binomial", "circus", "hendecahedron", "wide tent", "atipodal tent", "icosahedron binomial", "cube binomial", dodecahedron binomial
std::string policy;
    

void visualize_gauss_map(Forward3DSolver* forwardSolver){
  std::unique_ptr<ManifoldSurfaceMesh> sphere_mesh_ptr;
  std::unique_ptr<VertexPositionGeometry> sphere_geometry_ptr;
  std::tie(sphere_mesh_ptr, sphere_geometry_ptr) = generate_polyhedra("sphere");
  sphere_mesh = sphere_mesh_ptr.release();
  sphere_geometry = sphere_geometry_ptr.release();
  
  // update vis utils
  vis_utils.draw_gauss_map(forwardSolver, sphere_mesh, sphere_geometry);
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
    forwardSolver->hullGeometry->inputVertexPositions, forwardSolver->hullMesh->getFaceVertexList());
  forwardSolver->hullGeometry->requireFaceNormals();
  psHullMesh->addFaceVectorQuantity("face normals", forwardSolver->hullGeometry->faceNormals);
  psHullMesh->setEnabled(true);
  vis_utils.draw_G(forwardSolver);
}


void update_solver(){ // only doing this for convex input
  forwardSolver = new Forward3DSolver(mesh, geometry, G, true);
  forwardSolver->set_uniform_G();
  // TODO: doing this for hull dice search only?
  G = forwardSolver->get_G();
  forwardSolver->initialize_pre_computes();
  boundary_builder = new BoundaryBuilder(forwardSolver);
  boundary_builder->build_boundary_normals();
}


void generate_polyhedron_example(std::string poly_str, bool triangulate = false){
  // readManifoldSurfaceMesh()
  // std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(poly_str);
  std::tie(mesh_ptr, geometry_ptr) = generate_11_sided_polyhedron(poly_str);
  mesh = mesh_ptr.release();
  geometry = geometry_ptr.release();
  preprocess_mesh(mesh, geometry, triangulate, false, 1.);
}


void visualize_current_probs_and_goals(Forward3DSolver tmp_solver, 
                                       std::string policy_general, std::vector<std::pair<Vector3, double>> normal_prob_assignment, 
                                       bool show, bool print_probs = false){
  polyscope::registerPointCloud("Center of Mass", std::vector<Vector3>{tmp_solver.get_G()});
  auto curr_hull_psmesh = polyscope::registerSurfaceMesh("current hull", tmp_solver.hullGeometry->inputVertexPositions, tmp_solver.hullMesh->getFaceVertexList());
  
  curr_hull_psmesh->setSurfaceColor({0.1,0.9,0.1})->setEdgeWidth(2.)->setTransparency(0.7)->setEnabled(true);
  if (policy_general == "manual"){ // first word
    FaceData<double> my_probs = manual_stable_only_face_prob_assignment(&tmp_solver, normal_prob_assignment);
    curr_hull_psmesh->addFaceScalarQuantity("Goal probs", my_probs)->setColorMap("reds")->setEnabled(false);    
  }
  else if (policy_general == "manualCluster"){ // first word
    std::vector<std::tuple<std::vector<Face>, double, Vector3>> clustered_probs = manual_clustered_face_prob_assignment(&tmp_solver, normal_prob_assignment);
    FaceData<double> goal_cluster_probs(*tmp_solver.hullMesh, 0.),
                     current_cluster_probs(*tmp_solver.hullMesh, 0.);
    std::vector<Vector3> assignees;
    for (auto cluster: clustered_probs){
      double current_cluster_prob = 0.;
      std::vector<Face> faces = std::get<0>(cluster);
      double cluster_prob = std::get<1>(cluster);
      assignees.push_back(std::get<2>(cluster) + Vector3{0, 2, 0});
      for (Face f: faces){
        if (tmp_solver.face_last_face[f] == f){
          current_cluster_prob += tmp_solver.hullGeometry->faceArea(f)/(4.*PI);
          for (Face f2: tmp_solver.hullMesh->faces()){
            if (tmp_solver.face_last_face[f2] == f){
              goal_cluster_probs[f2] = cluster_prob;
            }
          }
        }
      }
      for (Face f: faces){
        current_cluster_probs[f] = current_cluster_prob;
      }
    }
    curr_hull_psmesh->addFaceScalarQuantity("Goal cluster probs", goal_cluster_probs)->setColorMap("reds")->setEnabled(false);    
    curr_hull_psmesh->addFaceScalarQuantity("current cluster accum probs", current_cluster_probs)->setColorMap("reds")->setEnabled(true); 
    polyscope::registerPointCloud("Cluster assignees", assignees);
  }

  BoundaryBuilder tmp_bnd_builder(&tmp_solver);
  tmp_bnd_builder.build_boundary_normals();
  vis_utils.update_visuals(&tmp_solver, &tmp_bnd_builder, sphere_mesh, sphere_geometry);
  FaceData<double> current_accum_probs(*tmp_solver.hullMesh, 0.);
  for (Face f: tmp_solver.hullMesh->faces()){
    current_accum_probs[f] = tmp_bnd_builder.face_region_area[tmp_solver.face_last_face[f]]/(4.*PI);
  }
  polyscope::getSurfaceMesh("current hull")->addFaceScalarQuantity("current probs", tmp_bnd_builder.face_region_area/(4.*PI))->setColorMap("reds")->setEnabled(false);    
  polyscope::getSurfaceMesh("current hull")->addFaceScalarQuantity("current accum probs", current_accum_probs)->setColorMap("reds")->setEnabled(true);    
  
  // Visuals/Logs
  if (print_probs){
    tmp_bnd_builder.print_area_of_boundary_loops();
  }
  if (show){
    // polyscope::frameTick();
    polyscope::screenshot(false);
    polyscope::show();
  }
}




void initialize_state(std::string input_name){
    generate_polyhedron_example(input_name);
    update_solver();
    init_visuals();
    visualize_gauss_map(forwardSolver);//
    boundary_builder->print_area_of_boundary_loops();
    vis_utils.update_visuals(forwardSolver, boundary_builder, sphere_mesh, sphere_geometry);
    // TODO : temporary
    Forward3DSolver tmp_solver(vertex_data_to_matrix(forwardSolver->hullGeometry->inputVertexPositions), 
                               vec32vec(forwardSolver->get_G()), true);
    tmp_solver.initialize_pre_computes();
    for (Face f: tmp_solver.hullMesh->faces()){
      std::cout << "face " << f.getIndex() << " N:" << tmp_solver.hullGeometry->faceNormal(f) << std::endl;
    }
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


void get_dice_energy_grads(Eigen::MatrixX3d hull_positions, Eigen::Vector3d G_vec, double bary_reg, double coplanar_reg, double cluster_distance_reg,
                           Eigen::MatrixX3d &df_dv, Eigen::Vector3d &df_dG, double &dice_energy,
                           bool use_autodiff, bool frozen_G, 
                           std::string policy_general, std::vector<std::pair<Vector3, double>> normal_prob_pairs, 
                           int fair_sides){
    
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
    // Eigen::MatrixX3d hull_positions = vertex_data_to_matrix(fwd_solver.hullGeometry->inputVertexPositions); 
    Forward3DSolver tmp_solver(hull_positions, G_vec, true); // indices shouldnt be shuffled here
    auto dice_energy_lambda = [&] <typename Scalar> (const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &hull_poses_G_append_vec) -> Scalar {
      // decompose flat vector to positions and center of mass; G is the last 3 elements
      Eigen::Vector3<Scalar> G_eigen = hull_poses_G_append_vec.tail(3);
      size_t flat_n = hull_poses_G_append_vec.rows();
      Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 3> > hull_poses(hull_poses_G_append_vec.head(flat_n-3).data(), flat_n/3 - 1, 3);
      return BoundaryBuilder::dice_energy<Scalar>(hull_poses, G_eigen, tmp_solver, 
                                                  bary_reg, coplanar_reg, cluster_distance_reg,
                                                  policy_general, normal_prob_pairs, fair_sides, true);
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


void dice_energy_opt(std::string policy, double bary_reg, double coplanar_reg, bool frozen_G, size_t step_count){
  polyscope::getSurfaceMesh("hull mesh")->setTransparency(0.5)->setEnabled(false);

  Forward3DSolver tmp_solver(mesh, geometry, G, true);
  tmp_solver.set_uniform_G();
  tmp_solver.initialize_pre_computes();
  
  G = tmp_solver.get_G();
  Eigen::Vector3d G_vec{G.x, G.y, G.z};
  Eigen::MatrixX3d hull_positions = vertex_data_to_matrix(tmp_solver.hullGeometry->inputVertexPositions);
  
  // inital assignment
  // policy shape
  std::string policy_general = policy.substr(0, policy.find(" ")); // first word
  std::string policy_shape = policy.substr(policy.find(" ") + 1); // second word
  std::cout << ANSI_FG_YELLOW << "policy general: " << policy_general << " policy shape: " << policy_shape << ANSI_RESET << "\n";
  std::vector<std::pair<Vector3, double>> normal_prob_pairs = normal_prob_assignment(policy_shape);
  
  // BoundaryBuilder tmp_bnd_builder(&tmp_solver);
  
  double current_sobolev_lambda = sobolev_lambda;
  double init_LS_step = dice_energy_step;
  double LS_step_tol = 1e-8;
  for (size_t iter = 0; iter < step_count; iter++){
    double dice_e;
    Eigen::Vector3d dfdG;
    Eigen::MatrixX3d dfdV(hull_positions.rows(), 3);
    
    printf("getting grads\n");
    FaceData<double> goal_probs;
    get_dice_energy_grads(hull_positions, G_vec, bary_reg, coplanar_reg, cluster_distance_reg,
                          dfdV, dfdG, dice_e, 
                          use_autodiff_for_dice_grad, frozen_G,
                          policy_general, normal_prob_pairs, fair_sides_count);
    // update_visuals_with_G(&tmp_solver, &tmp_bnd_builder);
    std::cout << ANSI_FG_YELLOW << "i: "<< iter << "\tDE: " << dice_e << ANSI_RESET << "\n";

    // DEBUG/visuals
    visualize_current_probs_and_goals(tmp_solver, policy_general, normal_prob_pairs, visualize_steps);
    // show grads
    polyscope::getSurfaceMesh("current hull")->addVertexVectorQuantity("ad total grads", -1.*dfdV)->setEnabled(true);
      
    // diffused grads
    if (do_sobolev_dice_grads){    
      // current_sobolev_lambda *= sobolev_lambda_decay;
      Eigen::MatrixXd diffused_dfdV = sobolev_diffuse_gradients(dfdV, *tmp_solver.hullMesh, current_sobolev_lambda, sobolev_p);
      
      double raw_norm = tinyAD_flatten(dfdV).norm(); // use eigen flaten
      double diffused_norm = tinyAD_flatten(diffused_dfdV).norm();
      
      // std::cout << "energies: " << energies[0][iter] << " " << energies[1][iter] << " " << energies[2][iter] << "\n";
      polyscope::getSurfaceMesh("current hull")->addVertexVectorQuantity(" sobolev diffused grads", diffused_dfdV*(-1.))->setEnabled(true);
      dfdV = diffused_dfdV;
    }

    // printf("line search\n");
    double opt_step_size = hull_update_line_search(dfdV, hull_positions, G_vec, bary_reg, coplanar_reg, cluster_distance_reg,
                                                   policy_general, normal_prob_pairs, fair_sides_count, 
                                                   init_LS_step, dice_search_decay, frozen_G, 1000, LS_step_tol);
    init_LS_step = opt_step_size; // TODO : adaptive step size

    std::cout << ANSI_FG_RED << "  line search step size: " << opt_step_size << ANSI_RESET << "\n";
    if (opt_step_size < LS_step_tol){
      if (!adaptive_reg){
        std::cout << ANSI_FG_RED << "  line search step size too small; breaking" << ANSI_RESET << "\n";
        break;
      }
      else{
        std::cout << ANSI_FG_RED << "  line search step size too small; modifying reg coeffs" << ANSI_RESET << "\n";
        std::cout << ANSI_FG_RED << "  bary reg: " << bary_reg << " coplanar reg: " << coplanar_reg << " sobolev lambda: " << current_sobolev_lambda << ANSI_RESET << "\n";
        init_LS_step = dice_energy_step; // reset if coeffs are changing
        // bary_reg /= dice_search_decay;
        // coplanar_reg /= dice_search_decay;

        current_sobolev_lambda *= sobolev_lambda_decay;
        if (bary_reg < 1e-3 || coplanar_reg < 1e-3 || bary_reg > 1e2 || coplanar_reg > 1e2 || current_sobolev_lambda < 0.1){
          std::cout << ANSI_FG_RED << "  regularizers too small/big; breaking" << ANSI_RESET << "\n";
          break;
        }
      }
    }
    hull_positions = hull_positions - opt_step_size * dfdV;

    tmp_solver = Forward3DSolver(hull_positions, G_vec, false); // could have been convaved
    if (!frozen_G){
      tmp_solver.set_uniform_G();
      G_vec = vec32vec(tmp_solver.get_G());
    }
    // IMPORTANT step; tmo solver's conv hull will shuffle the order of vertices
    hull_positions = vertex_data_to_matrix(tmp_solver.hullGeometry->inputVertexPositions);
    tmp_solver.initialize_pre_computes();
    if (policy_general != "fair"){
      std::vector<std::tuple<std::vector<Face>, double, Vector3>> clustered_face_normals = manual_clustered_face_prob_assignment(&tmp_solver, normal_prob_pairs);
      normal_prob_pairs = update_normal_prob_assignment(&tmp_solver, clustered_face_normals);
    }
  }

  // DEBUG/visuals
  visualize_current_probs_and_goals(tmp_solver, policy_general, normal_prob_pairs, false, true);

  optimized_mesh = tmp_solver.hullMesh;
  optimized_geometry = tmp_solver.hullGeometry; 

  // // DEBUG grad vecs
  // double dice_e;
  //        Eigen::Vector3d dfdG;
  //        Eigen::MatrixX3d dfdV(hull_positions.rows(), 3);
  // get_dice_energy_grads(hull_positions, G_vec, bary_reg, coplanar_reg,
  //                       dfdV, dfdG, dice_e, 
  //                       use_autodiff_for_dice_grad, frozen_G,
  //                       policy_general, normal_prob_pairs, fair_sides_count));
  // polyscope::getSurfaceMesh("current hull")->addVertexVectorQuantity("ad total grads", -1.*dfdV)->setEnabled(true);
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
  ImGui::SliderInt("ITERS",           &DE_step_count, 1, 200);
  ImGui::SliderInt("fair sides",      &fair_sides_count, 4, 20);
  ImGui::SliderFloat("DE step size",  &dice_energy_step, 0, 0.05);
  ImGui::SliderFloat("DE step decay", &dice_search_decay, 0.1, 1.);

  ImGui::SliderFloat("barycenter distance regularizer", &bary_reg, 0., 100.);
  ImGui::SliderFloat("coplanar regularizer", &coplanar_reg, 0., 100.);
  ImGui::SliderFloat("cluster distance regularizer", &cluster_distance_reg, 0., 100.);
  ImGui::Checkbox("sobolev grads", &do_sobolev_dice_grads);
  ImGui::SliderFloat("sobolev lambda", &sobolev_lambda, 0., 50.);
  ImGui::SliderFloat("decay sobolev lambda", &sobolev_lambda_decay, 0., 1.);
  ImGui::Checkbox("frozen G", &frozen_G);
  ImGui::Checkbox("adaptive reg", &adaptive_reg);
  ImGui::Checkbox("visualize steps", &visualize_steps);
  if (ImGui::Button("dice energy opt")){
    dice_energy_opt(policy, bary_reg, coplanar_reg, frozen_G, DE_step_count);
  }
  if (ImGui::Button("save optimized hull")){
    std::string output_name = input_name + "_optimized_dice";
    writeSurfaceMesh(*optimized_mesh, *optimized_geometry, "../meshes/hulls/" + std::string(output_name) +".obj");
  }
}


int main(int argc, char **argv) {
  // Parse args
  args::ArgumentParser parser("Dice energy optimization for convex shapes");
  // make parser for input name and other params
  args::HelpFlag help(parser,    "help", "Display this help menu", {'h', "help"});
  // args::Flag do_just_ours(parser, "do_just_ours", "do just ours", {"just_ours"});
  // args::ValueFlag<int> total_samples(parser, "ICOS_samples", "Total number of samples", {"samples"});
  args::ValueFlag<std::string> input_shape_arg(parser, "input_shape_str", "path to input shape", {'m', "mesh_dir"});
  args::ValueFlag<std::string> policy_arg(parser, "policy_general", " general policy string: fair | manual | manualCluster ", {'p', "policy"}, "manualCluster");


  try {
    parser.ParseCLI(argc, argv);
  }
  catch (args::Help) {
    std::cout << parser;
    return 0;
  }
  catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }
  catch (args::ValidationError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  if(input_shape_arg){
    input_name = args::get(input_shape_arg);
  }
  if(policy_arg){
    policy_general = args::get(policy_arg);
    // std::string policy_shape = "dodecahedron"; // "dodecahedron", "octahedron", "circus", "hendecahedron", "wide tent", "atipodal tent", "icosahedron", "cube"
    policy = policy_general + " " + input_name;
  }
  std::cout << ANSI_FG_YELLOW << "policy: " << policy << ANSI_RESET << "\n";
  std::cout << ANSI_FG_YELLOW << "input shape: " << input_name << ANSI_RESET << "\n";

  
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