#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "args.hxx"
#include "imgui.h"

#include "boundary_tools.h"
#include "inverse_design/dice_energy.h"
#include "mesh_factory.h"
// #include "geometry_utils.h" //; in fwd solver 
#include "visual_utils.h"
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

// DEBUG; will show every step in polyscope, go to next step by exiting polyscope window
bool visualize_steps = false; 


// probability assignment & policy
bool use_loaded_normal_probs = false;
std::string loaded_normal_prob_policy;
std::vector<std::pair<Vector3, double>> loaded_normal_prob_pairs;
std::vector<std::pair<Vector3, double>> normal_prob_pairs;
std::vector<std::pair<Vector3, double>> initial_normal_prob_pairs;

// log dice energy optimization steps
bool save_sequence_scr = false,
     save_sequence_files = false;
	 
std::string global_output_file = "../meshes/hulls/optimized.obj";



void load_mesh(const std::string& path, 
               ManifoldSurfaceMesh*& out_mesh, 
               VertexPositionGeometry*& out_geometry
			) {
    std::unique_ptr<SurfaceMesh> nm_mesh_ptr;
    std::unique_ptr<VertexPositionGeometry> nm_geometry_ptr;
    std::tie(nm_mesh_ptr, nm_geometry_ptr) = readSurfaceMesh(path);
    SurfaceMesh *nm_mesh = nm_mesh_ptr.release();
    VertexPositionGeometry *nm_geometry = nm_geometry_ptr.release();
    
    nm_mesh->greedilyOrientFaces();
    nm_mesh->compress();
    std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
    mesh_ptr = nm_mesh->toManifoldMesh();
    out_mesh = mesh_ptr.release();
    out_geometry = new VertexPositionGeometry(*out_mesh);
    
    // transfer from nm geometry
    for (Vertex v : out_mesh->vertices()) {
        out_geometry->inputVertexPositions[v.getIndex()] = nm_geometry->inputVertexPositions[v.getIndex()];
    }
    // Clean up temporary objects
    delete nm_mesh;
    delete nm_geometry;
}


void initialize_state(std::string input_path){
	load_mesh(input_path, mesh, geometry);
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
    // update_visuals(forwardSolver, boundary_builder, sphere_mesh, sphere_geometry);
}


void dice_energy_opt(
	std::string policy, double bary_reg, double coplanar_reg, 
	double cluster_distance_reg, double unstable_attraction_thresh,
	bool frozen_G, size_t step_count
){
	polyscope::getSurfaceMesh("init hull mesh")->setTransparency(0.5)->setEnabled(false);

	Forward3DSolver tmp_solver(mesh, geometry, G, true);
	tmp_solver.set_uniform_G();

	std::cout << ANSI_FG_YELLOW << "initializing for optimization" << ANSI_RESET<< "\n";
	tmp_solver.initialize_pre_computes();
	
	G = tmp_solver.get_G();
	Eigen::Vector3d G_vec{G.x, G.y, G.z};
	Eigen::MatrixX3d hull_positions = vertex_data_to_matrix(tmp_solver.hullGeometry->inputVertexPositions);
	
	// inital assignment
	// policy shape
	policy_general = policy.substr(0, policy.find(" ")); // first word
	std::string policy_shape = policy.substr(policy.find(" ") + 1); // second word
	std::cout << ANSI_FG_YELLOW << "policy general: " << policy_general << " policy shape: " << policy_shape << ANSI_RESET << "\n";

	if (use_loaded_normal_probs) {
		// use file-provided pairs and policy (if present)
		normal_prob_pairs = loaded_normal_prob_pairs;
		policy_general = loaded_normal_prob_policy;
		std::cout << ANSI_FG_GREEN << "Using loaded normal-probability assignment (policy: " << policy << ")" << ANSI_RESET << "\n";
	} else {
		// default behavior: generate assignment from built-in helpers
		normal_prob_pairs = normal_prob_assignment(policy_shape);
		if (policy_general == "fairCluster"){
			normal_prob_pairs = normal_prob_assignment_fair(&tmp_solver, dice_params.fair_sides_count);
			policy_general = "manualCluster"; // switch after initial assignment
		}
	}
	initial_normal_prob_pairs = normal_prob_pairs; // for saving to file later
	// BoundaryBuilder tmp_bnd_builder(&tmp_solver);
	
	double current_sobolev_lambda = dice_params.sobolev_lambda;
	double default_LS_step = dice_params.dice_energy_step;
	double step_tol = 1e-9;
	Eigen::MatrixX3d dfdV, diffused_dfdV;
	for (size_t iter = 0; iter < step_count; iter++){
		double dice_e;
		Eigen::Vector3d dfdG;
		dfdV = Eigen::MatrixX3d::Zero(hull_positions.rows(), 3);

		printf("getting grads\n");
		FaceData<double> goal_probs;
		get_dice_energy_grads(hull_positions, G_vec, bary_reg, coplanar_reg, cluster_distance_reg, unstable_attraction_thresh,
							dfdV, dfdG, dice_e, 
							frozen_G,
							policy_general, normal_prob_pairs, dice_params.fair_sides_count);
		std::cout << ANSI_FG_YELLOW << "i: "<< iter << "\tDE: " << dice_e << ANSI_RESET << "\n";

		// diffused grads
		if (dice_params.do_sobolev_dice_grads){
			current_sobolev_lambda *= dice_params.sobolev_lambda_decay;
			diffused_dfdV = sobolev_diffuse_gradients(dfdV, *tmp_solver.hullMesh, *tmp_solver.hullGeometry , 
														current_sobolev_lambda, dice_params.sobolev_p);      
			dfdV = diffused_dfdV;
		}
		// DEBUG/visuals
		std::cout << ANSI_FG_MAGENTA << "visualizing current probs and goals" << ANSI_RESET << "\n";
		visualize_current_probs_and_goals(tmp_solver, 
			sphere_mesh, sphere_geometry,
			policy_general, normal_prob_pairs, 
			dfdV, diffused_dfdV, 
			save_sequence_scr, save_sequence_files,
			visualize_steps, false, iter
		);
		printf("line search\n");
		double opt_step_size = 1.;
		if (dice_params.dice_search_decay < 1.){
			opt_step_size = hull_update_line_search(dfdV, hull_positions, G_vec, bary_reg, coplanar_reg, cluster_distance_reg, unstable_attraction_thresh,
															policy_general, normal_prob_pairs, dice_params.fair_sides_count, 
															default_LS_step, dice_params.dice_search_decay, frozen_G, 1000, step_tol);
		}

		std::cout << ANSI_FG_RED << "  line search step size: " << opt_step_size << ANSI_RESET << "\n";
		std::cout << ANSI_FG_MAGENTA << "  sobolev lambda: " << current_sobolev_lambda << ANSI_RESET << "\n";
		if (opt_step_size < step_tol){
			std::cout << ANSI_FG_RED << "  line search step size too small; breaking" << ANSI_RESET << "\n";
			break;
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
		if (policy_general != "fair"){ // dynamic assignment; changes at each iteration
			std::vector<std::tuple<std::vector<Face>, double, Vector3>> clustered_face_normals = manual_clustered_face_prob_assignment(&tmp_solver, normal_prob_pairs);
			normal_prob_pairs = update_normal_prob_assignment(&tmp_solver, clustered_face_normals, dice_params.update_with_max_prob_face);
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
	if (ImGui::SliderInt("fair sides",      &dice_params.fair_sides_count, 2, 20));
		// fair_sides_count = dice_params.fair_sides_count;
	if (ImGui::SliderFloat("DE step size",  &dice_params.dice_energy_step, 0, 0.05));
		// dice_energy_step = dice_params.dice_energy_step;
	if(ImGui::SliderFloat("DE line seach decay", &dice_params.dice_search_decay, 0.1, 1.));
		// dice_search_decay = dice_params.dice_search_decay;
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
			dice_params.cluster_distance_reg, dice_params.unstable_attraction_thresh,
						dice_params.frozen_G, dice_params.DE_step_count);
	}
	if (ImGui::Button("save optimized hull")){
		// Ensure output directory exists
		ensure_directory_exists(global_output_file);
		
		// Save optimized mesh
		writeSurfaceMesh(*optimized_mesh, *optimized_geometry, global_output_file);
		std::cout << "Optimized mesh saved to: " << global_output_file << std::endl;
		
		// Generate names for parameter and G files in the same directory
		std::filesystem::path out_path(global_output_file);
		std::string base_name = out_path.stem().string(); // filename without extension
		std::string dir_path = out_path.parent_path().string();
		
		std::string params_file = dir_path + "/params.json";
		std::string G_file_path = dir_path + "/G_" + base_name + ".txt";
		
		// // Save initial mesh
		// std::string initial_mesh_file = dir_path + "/initial_mesh.obj";
		// writeSurfaceMesh(*mesh, *geometry, initial_mesh_file);
		// std::cout << "Initial mesh saved to: " << initial_mesh_file << std::endl;
		
		// Save parameters
		save_convex_dice_params(dice_params, params_file);
		
		if (dice_params.frozen_G) { // otherwise uniform mass is assumed
			// Save G to separate text file (for backward compatibility)
			std::ofstream G_file(G_file_path);
			G_file << G.x << " " << G.y << " " << G.z << "\n";
			G_file.close();
			std::cout << "Center of mass saved to: " << G_file_path << std::endl;
		}

		// save normal-probability pairs
		std::string np_file_path = dir_path + "/final_normal_probs.json";
		save_vector3_double_pairs(normal_prob_pairs, np_file_path, policy_general);
		std::cout << "Normal-probability pairs saved to: " << np_file_path << "\n";
		// save initial normal-probability pairs
		std::string initial_np_file_path = dir_path + "/initial_normal_probs.json";
		save_vector3_double_pairs(initial_normal_prob_pairs, initial_np_file_path, policy_general);
		std::cout << "Initial normal-probability pairs saved to: " << initial_np_file_path << "\n";
	}
}


int main(int argc, char **argv) {
	// Parse args
	args::ArgumentParser parser("Dice energy optimization for convex shapes");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	args::ValueFlag<std::string> input_shape_arg(parser, "input_shape_str", "path to input shape", {'m', "mesh"});
	args::ValueFlag<std::string> policy_arg(parser, "policy_general", "general policy string: fair | manual | manualCluster | fairCluster", {'p', "policy"}, "manualCluster");
	args::ValueFlag<std::string> param_file_arg(parser, "param_file", "path to parameter JSON file to load", {"params"});
	args::ValueFlag<std::string> out_file_arg(parser, "output_file", "path to output optimized mesh (.obj)", {'o', "out"});
	args::ValueFlag<std::string> normal_prob_file_arg(parser, "normal_prob_file", "path to normal-probability JSON file (contains pairs + policy)", {"normal_probs"});


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

	// Default output file if not provided
	std::string output_file = "../meshes/hulls/optimized.obj";
	if (out_file_arg) {
		output_file = args::get(out_file_arg);
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
	// After parsing and after loading dice params, load normal probs if requested
  	if (normal_prob_file_arg) {
		std::string np_path = args::get(normal_prob_file_arg);
		std::string tmp_policy;
		auto pairs = load_vector3_double_pairs(np_path, &tmp_policy);
		if (!pairs.empty()) {
			loaded_normal_prob_pairs = std::move(pairs);
			loaded_normal_prob_policy = tmp_policy;
			use_loaded_normal_probs = true;
			std::cout << "Loaded normal-probability assignment from: " << np_path << "\n";
			if (!loaded_normal_prob_policy.empty())
			std::cout << "Loaded policy: " << loaded_normal_prob_policy << "\n";
		} else {
			std::cerr << "Warning: no pairs loaded from " << np_path << " — continuing without it\n";
		}
	}

	std::cout << ANSI_FG_YELLOW << "policy: " << policy << ANSI_RESET << "\n";
	std::cout << ANSI_FG_YELLOW << "input shape: " << input_name << ANSI_RESET << "\n";
	std::cout << ANSI_FG_YELLOW << "output file: " << output_file << ANSI_RESET << "\n";
	global_output_file = output_file;

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