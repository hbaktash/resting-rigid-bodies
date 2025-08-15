#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "args.hxx"
#include "imgui.h"

#include "boundary_tools.h"
#include "mesh_factory.h"
// #include "geometry_utils.h" //; in fwd solver 
#include "visual_utils.h"
#include <stan/math.hpp>
#include <nlohmann/json.hpp>

#include "file_IO.h"

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
std::string input_name;
std::string policy_general;
std::string policy;
Vector3 G = Vector3::zero();  

VisualUtils vis_utils;

ManifoldSurfaceMesh* sphere_mesh;
VertexPositionGeometry* sphere_geometry;


// Replace all global parameter variables with a single struct
ConvexDiceParams dice_params;


int fair_sides_count = 6, // for optimization
    DE_step_count = 40;
bool do_sobolev_dice_grads = false,
     frozen_G = false,
     use_autodiff_for_dice_grad = true,
     visualize_steps = true,
     adaptive_reg = false;
float sobolev_lambda = 2.,
      sobolev_lambda_decay = 0.8,
      dice_energy_step = 0.01,
      dice_search_decay = 0.98,
      bary_reg = 0.1,
      coplanar_reg = 0.0,
      cluster_distance_reg = 0.0,
      unstable_attraction_thresh = 0.1;
int sobolev_p = 2;
// optimization stuff
bool update_with_max_prob_face = true;

// log stuff
bool save_sequence_scr = false,
     save_sequence_files = false;

// // example choice
// std::vector<std::string> all_input_names = {std::string("6 prism"), std::string("hendecahedron"), std::string("triangular"), std::string("circus"), std::string("icosahedron"), std::string("dodecahedron"), std::string("cuub"), std::string("octahedron")}; // {std::string("tet"), std::string("tet2"), std::string("cube"), std::string("tilted cube"), std::string("dodecahedron"), std::string("Conway spiral 4"), std::string("oloid")};



void generate_polyhedron_example(std::string poly_str){
  std::tie(mesh_ptr, geometry_ptr) = generate_11_sided_polyhedron(poly_str);
  mesh = mesh_ptr.release();
  geometry = geometry_ptr.release();
}


void initialize_state(std::string input_name){
	std::cout << "loaded input name: " << input_name << std::endl;
    generate_polyhedron_example(input_name);
    bool triangulate = false;
    preprocess_mesh(mesh, geometry, triangulate, false);
    center_and_normalize(mesh, geometry);
    // set global G to uniform here
	Forward3DSolver* forwardSolver = new Forward3DSolver(mesh, geometry, G, true);
    forwardSolver->set_uniform_G();
    G = forwardSolver->get_G();
    std::cout << ANSI_FG_GREEN << "G is set to uniform:" << G << ANSI_RESET << std::endl;
    forwardSolver->initialize_pre_computes();
    BoundaryBuilder *boundary_builder = new BoundaryBuilder(forwardSolver);
    boundary_builder->build_boundary_normals();

    // visuals
    init_visuals(mesh, geometry, forwardSolver, boundary_builder);
    update_visuals(forwardSolver, boundary_builder, sphere_mesh, sphere_geometry);
}


void get_dice_energy_grads(Eigen::MatrixX3d hull_positions, Eigen::Vector3d G_vec, double bary_reg, double coplanar_reg, double cluster_distance_reg, double unstable_attraction_thresh,
                           Eigen::MatrixX3d &df_dv, Eigen::Vector3d &df_dG, double &dice_energy,
                           bool use_autodiff, bool frozen_G, 
                           std::string policy_general, std::vector<std::pair<Vector3, double>> normal_prob_pairs, 
                           int fair_sides){
  // Eigen::MatrixX3d hull_positions = vertex_data_to_matrix(fwd_solver.hullGeometry->inputVertexPositions); 
  Forward3DSolver tmp_solver(hull_positions, G_vec, true); // indices shouldnt be shuffled here
  auto dice_energy_lambda = [&] <typename Scalar> (const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &hull_poses_G_append_vec) -> Scalar {
    // decompose flat vector to positions and center of mass; G is the last 3 elements
    Eigen::Vector3<Scalar> G_eigen = hull_poses_G_append_vec.tail(3);
    size_t flat_n = hull_poses_G_append_vec.rows();
    Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 3> > hull_poses(hull_poses_G_append_vec.head(flat_n-3).data(), flat_n/3 - 1, 3);
    return BoundaryBuilder::dice_energy<Scalar>(hull_poses, G_eigen, tmp_solver, 
                                                bary_reg, coplanar_reg, cluster_distance_reg, unstable_attraction_thresh,
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


void dice_energy_opt(std::string policy, double bary_reg, double coplanar_reg, bool frozen_G, size_t step_count){
  polyscope::getSurfaceMesh("init hull mesh")->setTransparency(0.5)->setEnabled(false);

  Forward3DSolver tmp_solver(mesh, geometry, G, true);
  tmp_solver.set_uniform_G();

  std::cout << ANSI_FG_YELLOW << "initializing for optimization" << ANSI_RESET<< "\n";
  // std::cout << " G was: "<< G << "\n";
  // std::cout << " G is: "<< tmp_solver.get_G() << "\n";
  // polyscope::show();
  tmp_solver.initialize_pre_computes();
  
  G = tmp_solver.get_G();
  Eigen::Vector3d G_vec{G.x, G.y, G.z};
  Eigen::MatrixX3d hull_positions = vertex_data_to_matrix(tmp_solver.hullGeometry->inputVertexPositions);
  
  // inital assignment
  // policy shape
  std::string policy_general = policy.substr(0, policy.find(" ")); // first word
  std::string policy_shape = policy.substr(policy.find(" ") + 1); // second word
  std::cout << ANSI_FG_YELLOW << "policy general: " << policy_general << " policy shape: " << policy_shape << ANSI_RESET << "\n";
  // if (policy_shape == "circus" || policy_shape == "hendecahedron" || policy_shape == "cube"|| policy_shape == "dodecahedron"|| policy_shape == "octahedron"|| policy_shape.substr(policy_shape.find(" ") + 1) == "prism")
  std::vector<std::pair<Vector3, double>> normal_prob_pairs;
  normal_prob_pairs = normal_prob_assignment(policy_shape);
  if (policy_general == "fairCluster"){
    normal_prob_pairs = normal_prob_assignment_fair(&tmp_solver, fair_sides_count);
    policy_general = "manualCluster"; // switch after initial assignment
  }
  // BoundaryBuilder tmp_bnd_builder(&tmp_solver);
  
  double current_sobolev_lambda = sobolev_lambda;
  double init_LS_step = dice_energy_step;
  double LS_step_tol = 1e-9;
  Eigen::MatrixX3d dfdV, diffused_dfdV;
  for (size_t iter = 0; iter < step_count; iter++){
    double dice_e;
    Eigen::Vector3d dfdG;
    dfdV = Eigen::MatrixX3d::Zero(hull_positions.rows(), 3);

    printf("getting grads\n");
    FaceData<double> goal_probs;
    get_dice_energy_grads(hull_positions, G_vec, bary_reg, coplanar_reg, cluster_distance_reg, unstable_attraction_thresh,
                          dfdV, dfdG, dice_e, 
                          use_autodiff_for_dice_grad, frozen_G,
                          policy_general, normal_prob_pairs, fair_sides_count);
    // update_visuals_with_G(&tmp_solver, &tmp_bnd_builder);
    std::cout << ANSI_FG_YELLOW << "i: "<< iter << "\tDE: " << dice_e << ANSI_RESET << "\n";

    // diffused grads
    if (do_sobolev_dice_grads){
      current_sobolev_lambda *= sobolev_lambda_decay;
      diffused_dfdV = sobolev_diffuse_gradients(dfdV, *tmp_solver.hullMesh, *tmp_solver.hullGeometry , 
                                                current_sobolev_lambda, sobolev_p);      
      dfdV = diffused_dfdV;
    }
    // DEBUG/visuals
    visualize_current_probs_and_goals(tmp_solver, 
      sphere_mesh, sphere_geometry,
      policy_general, normal_prob_pairs, 
      dfdV, diffused_dfdV, 
      save_sequence_scr, save_sequence_files,
      visualize_steps, false, iter);
    // printf("line search\n");
    double opt_step_size = 1.;
    if (dice_search_decay != 1.){
      opt_step_size = hull_update_line_search(dfdV, hull_positions, G_vec, bary_reg, coplanar_reg, cluster_distance_reg, unstable_attraction_thresh,
                                                    policy_general, normal_prob_pairs, fair_sides_count, 
                                                    init_LS_step, dice_search_decay, frozen_G, 1000, LS_step_tol);
      init_LS_step = opt_step_size < dice_energy_step/100. ? opt_step_size * 20 : opt_step_size; // TODO : adaptive step size
    }

    std::cout << ANSI_FG_RED << "  line search step size: " << opt_step_size << ANSI_RESET << "\n";
    std::cout << ANSI_FG_MAGENTA << "  sobolev lambda: " << current_sobolev_lambda << ANSI_RESET << "\n";
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
      normal_prob_pairs = update_normal_prob_assignment(&tmp_solver, clustered_face_normals, update_with_max_prob_face);
    }
  }

  // DEBUG/visuals
  visualize_current_probs_and_goals(
    tmp_solver, 
    sphere_mesh, sphere_geometry,
    policy_general, normal_prob_pairs, dfdV, diffused_dfdV, 
    save_sequence_scr, save_sequence_files,
    false, true, step_count);

  optimized_mesh = tmp_solver.hullMesh;
  optimized_geometry = tmp_solver.hullGeometry; 
}


// Update myCallback() to use the struct
void myCallback() {
  ImGui::SliderInt("ITERS",           &dice_params.DE_step_count, 1, 200);
  ImGui::SliderInt("fair sides",      &dice_params.fair_sides_count, 2, 20);
  ImGui::SliderFloat("DE step size",  &dice_params.dice_energy_step, 0, 0.05);
  ImGui::SliderFloat("DE step decay", &dice_params.dice_search_decay, 0.1, 1.);

  ImGui::SliderFloat("barycenter distance regularizer", &dice_params.bary_reg, 0., 100.);
  ImGui::SliderFloat("coplanar regularizer", &dice_params.coplanar_reg, 0., 100.);
  ImGui::SliderFloat("cluster distance regularizer", &dice_params.cluster_distance_reg, 0., 100.);
  ImGui::SliderFloat("cluster unstable-stable attraction threshold", &dice_params.unstable_attraction_thresh, 0., 100.);

  ImGui::Checkbox("sobolev grads", &dice_params.do_sobolev_dice_grads);
  ImGui::SliderFloat("sobolev lambda", &dice_params.sobolev_lambda, 0., 50.);
  ImGui::SliderFloat("decay sobolev lambda", &dice_params.sobolev_lambda_decay, 0., 1.);
  ImGui::Checkbox("frozen G", &dice_params.frozen_G);
  ImGui::Checkbox("update clusterN with max prob face", &dice_params.update_with_max_prob_face);
  ImGui::Checkbox("adaptive reg", &dice_params.adaptive_reg);
  ImGui::Checkbox("visualize steps", &visualize_steps);
  ImGui::Checkbox("save scrs", &save_sequence_scr);
  ImGui::Checkbox("save files", &save_sequence_files);

  if (ImGui::Button("dice energy opt")){
    dice_energy_opt(policy, dice_params.bary_reg, dice_params.coplanar_reg, 
                    dice_params.frozen_G, dice_params.DE_step_count);
  }
  if (ImGui::Button("save optimized hull")){
    std::string output_name = input_name + "_d" + std::to_string(dice_params.fair_sides_count) + "_" + policy_general;
    writeSurfaceMesh(*optimized_mesh, *optimized_geometry, "../meshes/hulls/" + output_name + ".obj");
    
    // Update struct with current G and save parameters
    save_convex_dice_params(dice_params, "../meshes/hulls/params_" + output_name + ".json");

    // save G to separate text file (for backward compatibility)
    std::ofstream G_file;
    G_file.open("../meshes/hulls/G_" + output_name + ".txt");
    G_file << G.x << " " << G.y << " " << G.z << "\n";
    G_file.close();
  }
}

int main(int argc, char **argv) {
  // Parse args
  args::ArgumentParser parser("Dice energy optimization for convex shapes");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::ValueFlag<std::string> input_shape_arg(parser, "input_shape_str", "path to input shape", {'m', "mesh_dir"});
  args::ValueFlag<std::string> policy_arg(parser, "policy_general", "general policy string: fair | manual | manualCluster | fairCluster", {'p', "policy"}, "manualCluster");
  args::ValueFlag<std::string> param_file_arg(parser, "param_file", "path to parameter JSON file to load", {"params"});

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

  // Load parameters from file if provided
  if (param_file_arg) {
    if (!load_convex_dice_params(dice_params, args::get(param_file_arg))) {
      std::cerr << "Failed to load parameters, using defaults" << std::endl;
    }
  }

  if (input_shape_arg) {
    input_name = args::get(input_shape_arg);
  }
  if (policy_arg) {
    policy_general = args::get(policy_arg);
    policy = policy_general + " " + input_name;
  }

  std::cout << ANSI_FG_YELLOW << "policy: " << policy << ANSI_RESET << "\n";
  std::cout << ANSI_FG_YELLOW << "input shape: " << input_name << ANSI_RESET << "\n";


  // build mesh
  vis_utils = VisualUtils();

  polyscope::init();
  initialize_state(input_name);
  input_name = input_name.substr(input_name.find_last_of("/") + 1);
  input_name = input_name.substr(0, input_name.find_last_of("."));
  
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