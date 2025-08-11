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
// #define BT_USE_DOUBLE_PRECISION

#include "unsupported/Eigen/EulerAngles"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
// #include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include "nlohmann/json.hpp"
// #include "polyscope/nlohmann/json.hpp"
// #include "bullet3/examples/BasicExample.h"
#include "args.hxx"
#include "imgui.h"

#include "coloring.h"
#include "forward3D.h"
#include "mesh_factory.h"
#include "geometry_utils.h"
// #include "bullet_sim.h"
#include "visual_utils.h"

// bullet stuff
#include "LinearMath/btVector3.h"
#include "btBulletDynamicsCommon.h"
#include "bullet_sim.h"

// #include "ipc/ipc.hpp"

// system stuff
#include "chrono"
#include <filesystem>
#include <fstream> 

namespace fs = std::filesystem;

namespace chrono = std::chrono;
using clock_type = chrono::high_resolution_clock;
using seconds_fp = chrono::duration<double, chrono::seconds::period>;

using namespace geometrycentral;
using namespace geometrycentral::surface;

// simulation stuff
// Bullet simulation stuff
PhysicsEnv* my_env;

Vector3 G;

float orientation_gui[3] = {0., -1., 0.};


// example choice (removed unused all_polyhedra_items/all_polygons_current_item)

// GC stuff
std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr;
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

bool verbose = false,
     bullet_sim = false,
     ipc_sim = false,
     just_ours = false;

// raster image stuff
FaceData<Vector3> face_colors;
int ICOS_samples = 10, 
    max_steps_IPC = 2000;
bool ICOS_sampling = true;

// quasi static simulation stuff
Forward3DSolver* forwardSolver;
BoundaryBuilder *boundary_builder;
VisualUtils vis_utils;


int snail_trail_dummy_counter = 0;
int solid_angle_face_index = 0;


// for file saving
std::string mesh_title = "mesh";
bool save_oriented = false,
     save_orientation_trail_to_file = false,
     save_orientation_trail_scrs = false,
     snapshot_trail = false,
     save_quasi_trail_to_file = false,
     snapshot_quasi_trail = false;

int stable_face_index = 0; // index of the stable face to visualize
bool save_stable_orientation_pose = false;

std::string multi_orientations_path;
// Functions 


void polyscope_defaults(){
	// hide everything
    polyscope::getSurfaceMesh("init input mesh")->setEnabled(false);
    polyscope::getSurfaceMesh("init hull mesh")->setEnabled(false);
    polyscope::getCurveNetwork("Arc curves all edge arcs")->setEnabled(false);
    polyscope::getCurveNetwork("Arc curves region boundaries")->setEnabled(false);
    polyscope::getPointCloud("Edge equilibria")->setEnabled(false);
    polyscope::getPointCloud("stable Face Normals")->setEnabled(false);
    polyscope::getPointCloud("stable Vertices Normals")->setEnabled(false);
    polyscope::getPointCloud("Center of Mass")->setEnabled(false);

    polyscope::options::ssaaFactor = 3; // supersampling
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
}


void load_mesh(std::string path, bool preprocess = true){
	std::unique_ptr<SurfaceMesh> nm_mesh_ptr;
	std::unique_ptr<VertexPositionGeometry> nm_geometry_ptr;
	std::tie(nm_mesh_ptr, nm_geometry_ptr) = readSurfaceMesh(path);
	SurfaceMesh *nm_mesh = nm_mesh_ptr.release();
	VertexPositionGeometry *nm_geometry = nm_geometry_ptr.release();
	
	nm_mesh->greedilyOrientFaces();
	nm_mesh->compress();
	mesh_ptr = nm_mesh->toManifoldMesh();
	mesh = mesh_ptr.release();
	geometry = new VertexPositionGeometry(*mesh);
	// transfer from nm geometry
	for (Vertex v : mesh->vertices()) {
		geometry->inputVertexPositions[v.getIndex()] = nm_geometry->inputVertexPositions[v.getIndex()];
	}

	// preproccess and shift for external use
	if (preprocess){
		preprocess_mesh(mesh, geometry, true, false);
		G = find_center_of_mass(*mesh, *geometry).first;
		// std::cout << "center of mass before shift: " << G << "\n";
		for (Vertex v: mesh->vertices()){
		    geometry->inputVertexPositions[v] -= G;
		}
		G = find_center_of_mass(*mesh, *geometry).first;
	}
	G = find_center_of_mass(*mesh, *geometry).first;
}


void update_solver(){
  forwardSolver = new Forward3DSolver(mesh, geometry, G, true);
  forwardSolver->initialize_pre_computes();
  boundary_builder = new BoundaryBuilder(forwardSolver);
  boundary_builder->build_boundary_normals();
}


void initialize_bullet_env(Vector3 orientation, bool visuals = true){
    // params 
    double bullet_step_size = 0.01; // 0.016 is the default for 60 FPS;
    double ground_box_y = -2.1;
    Vector3 ground_box_shape({10,1,10});
    // physics env
    my_env = new PhysicsEnv();
    my_env->init_physics();
    my_env->init_geometry(forwardSolver->hullMesh, forwardSolver->hullGeometry);
    my_env->add_ground(ground_box_y, ground_box_shape);
    my_env->add_object(G, orientation);
    my_env->default_step_size = bullet_step_size;
}

std::pair<std::vector<Eigen::Matrix4d>, std::vector<Vector3>> 
generate_transformations_for_orientation_sequence(Vector3 initial_orientation,
    Forward3DSolver* forwardSolver, Vector3 floor_vec, double goal_angle_step = 0.005){
    // Get the snail trail and make forward log
    std::vector<Vector3> snail_trail = forwardSolver->snail_trail_log(initial_orientation);
    // split the snail trail
    std::vector<Vector3> snail_trail_refined;
    for (int i = 1; i < snail_trail.size()-1; i++){
        Vector3 local_axis = cross(snail_trail[i-1], snail_trail[i]).normalize();
        double local_total_angle = angle(snail_trail[i-1], snail_trail[i]);
        int steps = (int)ceil(local_total_angle/goal_angle_step) + 1;
        // in steps = 2;
        for (int t = 0; t < steps; t++){
            double angle_0 = local_total_angle * (double)t/double(steps);
            Vector3 normal_0 = snail_trail[i-1].rotateAround(local_axis, angle_0);
            snail_trail_refined.push_back(normal_0);
        }
    }
    snail_trail_refined.push_back(snail_trail[snail_trail.size()-1]);

    VertexData<Vector3> init_hull_positions = forwardSolver->hullGeometry->inputVertexPositions;
    VertexData<Vector3> tmp_hull_positions(*forwardSolver->hullMesh);
    tmp_hull_positions = init_hull_positions;
    
    std::vector<Vector3> saved_snail_trail_refined = snail_trail_refined; // for future output
    // get transformations
    std::vector<Eigen::Matrix4d> transformations;
    // transformations.push_back(Eigen::Matrix4d::Identity()); // initial identity matrix
    for (int i = 0; i < snail_trail_refined.size(); i++){
        // get n_i locally
        // SO3 conversion
        Vector3 normal_0 = snail_trail_refined[i];
        Vector3 rot_axis = cross(normal_0, floor_vec).normalize();
        double rot_angle = angle(normal_0, floor_vec);

        for (int j = 0; j < snail_trail_refined.size(); j++){
            snail_trail_refined[j] = snail_trail_refined[j].rotateAround(rot_axis, rot_angle);
        }
        // shift contact to origin
        double lowest_height = 1e4;
        Vector3 contact_p;
        for (Vertex v: forwardSolver->hullMesh->vertices()){
            if (tmp_hull_positions[v].y < lowest_height){
                lowest_height = tmp_hull_positions[v].y;
                contact_p = tmp_hull_positions[v];
            }
        }
        Eigen::AngleAxisd aa(rot_angle, Eigen::Vector3d(rot_axis.x, rot_axis.y, rot_axis.z));
        Eigen::Matrix3d rotation_matrix = aa.toRotationMatrix();
        // do the rotation around the contact point
        // shifting contact
        tmp_hull_positions -= contact_p;
        for (Vertex v: forwardSolver->hullMesh->vertices()){
            // tmp_hull_positions[v] = tmp_hull_positions[v].rotateAround(rot_axis, rot_angle);
            Eigen::Vector3d tmp_v(tmp_hull_positions[v].x, tmp_hull_positions[v].y, tmp_hull_positions[v].z);
            tmp_v = rotation_matrix * tmp_v;
            tmp_hull_positions[v] = Vector3({tmp_v[0], tmp_v[1], tmp_v[2]});
        }
        // shift back
        tmp_hull_positions += contact_p;
        // axis angle rotation to matrix

        // correct height; not really necessary for most steps; only when the contact point is a triangle, and the wrong vertex is chosen, can be fixed
        lowest_height = 1e4;
        for (Vertex v: forwardSolver->hullMesh->vertices()){
            if (tmp_hull_positions[v].y < lowest_height)
                lowest_height = tmp_hull_positions[v].y;
        }
        Eigen::Matrix4d transformation_matrix;
        transformation_matrix.setIdentity();
        transformation_matrix.block<3,3>(0,0) = rotation_matrix;
        // the last column is the translation
        Eigen::Vector3d contact_p_eigen(contact_p.x, contact_p.y, contact_p.z);
        Eigen::Vector3d translation = -rotation_matrix * contact_p_eigen + contact_p_eigen + Eigen::Vector3d(0, -lowest_height, 0);
        transformation_matrix.block<3,1>(0,3) = translation;
        transformations.push_back(transformation_matrix);
        // update tmp_hull_positions height
        for (Vertex v: forwardSolver->hullMesh->vertices()){
            tmp_hull_positions[v].y += -lowest_height;
        }
    }
    return {transformations, saved_snail_trail_refined};
}


std::pair<std::vector<Eigen::Matrix4d> , std::vector<Vector3>>
quasi_static_snail_trail_to_global_transformations(Forward3DSolver* forwardSolver, Vector3 initial_orientation, double refinement_step_size, bool visualize = false){
    // initialize the trail in forwardSolver
    std::vector<Vector3> snail_trail = forwardSolver->snail_trail_log(initial_orientation);
    // split the snail trail
    auto [local_transformation_matrices, saved_snail_trail_refined] = generate_transformations_for_orientation_sequence(
        initial_orientation, forwardSolver, Vector3({0,-1,0}), refinement_step_size);
    // local to global; L_n = R_n * R_(n-1) * ... * R_0
    std::vector<Eigen::Matrix4d> transformation_matrices;
    Eigen::Matrix4d global_transformation = Eigen::Matrix4d::Identity();
    for (int i = 0; i < local_transformation_matrices.size(); i++){
        Eigen::Matrix4d local_transformation = local_transformation_matrices[i];
        global_transformation = local_transformation * global_transformation;
        transformation_matrices.push_back(global_transformation);
    }
    // apply to the mesh and visualize each step
    if (visualize){
        visualize_quasi_static_drop_sequence(transformation_matrices, saved_snail_trail_refined, forwardSolver);
    }
    return {transformation_matrices, saved_snail_trail_refined};
}


VertexData<Vector3> orient_mesh_with_down_vec(Vector3 orientation, Vector3 down_vec,
                               ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry){
    Vector3 axis = cross(orientation, down_vec);
    if (axis.norm() < 1e-6) // parallel
        axis = Vector3({1,0,0});
    axis = axis.normalize();
    double angle = acos(dot(orientation, down_vec));
    VertexData<Vector3> new_positions(*mesh);
    for (Vertex v: mesh->vertices()){
        Vector3 p = geometry->inputVertexPositions[v];
        new_positions[v] = p.rotateAround(axis, angle);
    }
    return new_positions;
}


// polyscope callback
void myCallback() {
    // sliders for orientation
    ImGui::SliderFloat3("orientation", orientation_gui, -10.0f, 10.0f);
    if (ImGui::Button("show orientation and QS trail")) {
        // use-case 1; show orientation as a point on S^2 and orient the object accordingly, later save the oriented obj to file
        Vector3 orientation({orientation_gui[0], orientation_gui[1], orientation_gui[2]});
        orientation = orientation.normalize();
        // Make a rotated mesh from the rotation that rotates {0,-1,0} to orientation
        Vector3 down_vec({0,-1,0});
        VertexData<Vector3> new_positions = orient_mesh_with_down_vec(orientation, down_vec, mesh, geometry);

        polyscope::registerSurfaceMesh("oriented mesh", new_positions, mesh->getFaceVertexList());
        // show the orientation on the gauss map
        std::vector<Vector3> orientation_vec = {orientation * vis_utils.gm_radi + vis_utils.center};
        // black rgb
        polyscope::registerPointCloud("orientation on gm", orientation_vec)->setPointRadius(0.03, false)->setPointColor({0,0,0});
    }
    if (ImGui::Button("Build quasistatic snail trail")){
        // get rotations corresponding to the snail trail
        Vector3 orientation({orientation_gui[0], orientation_gui[1], orientation_gui[2]});
        orientation = orientation.normalize();
        
        // initialize bullet env
        std::vector<Vector3> snail_trail = forwardSolver->snail_trail_log(orientation);
        draw_trail_on_gm(snail_trail, {39./255., 189./255., 0}, "quasi-static trail", 3.);
        double refinement_step_size = 0.005;
        auto [transformation_matrices, saved_snail_trail_refined] = 
                quasi_static_snail_trail_to_global_transformations(
                    forwardSolver, orientation, refinement_step_size, false
                );
        
    }
    ImGui::Checkbox("save orientation trail", &save_orientation_trail_to_file);
    if (ImGui::Button("bullet snail trail")) {
        // initialize bullet env
        Vector3 orientation({orientation_gui[0], orientation_gui[1], orientation_gui[2]});
        orientation = orientation.normalize();
        initialize_bullet_env(orientation, false);
        // visualize snail trail 
        Face touching_face = my_env->final_stable_face(true);
        // animate the shape drop with Bullet
        std::vector<geometrycentral::DenseMatrix<double>> trans_mat_trail = my_env->trans_mat_trail;
        Vector<Vector3> init_positions = forwardSolver->inputGeometry->inputVertexPositions.toVector();
        // // save positions of the shape at every transformation
        std::vector<Vector<Vector3>> pos_trail;
        for (geometrycentral::DenseMatrix<double> tmp_trans: trans_mat_trail) {
            Vector<Vector3> tmp_positions = apply_trans_to_positions(init_positions, tmp_trans);
            pos_trail.push_back(tmp_positions);
        }
        
        // draw snail trail on gm
        std::vector<Vector3> snail_trail = my_env->orientation_trail;
        for (Vector3 &v: snail_trail)
            v += vis_utils.center;
        auto trail_pc = polyscope::registerPointCloud("bullet pc trail", snail_trail);
        glm::vec3 init_color = {0.8,0.8,0.2};
        // RGB is 209 227 28
        glm::vec3 bullet_color = {209./255., 227./255., 28./255.};
        trail_pc->setPointColor(bullet_color);
        trail_pc->setPointRadius(0.003, false);
        // trail_pc->setEnabled(true);
        // show the last position of the snail trail with red
        std::vector<Vector3> last_snail_trail = {snail_trail.back()};
        polyscope::registerPointCloud("final bullet trail", last_snail_trail)->setPointColor({1,0,0})->setPointRadius(0.03, false);

        // height function at every orientation
        std::vector<double> height_log;
        for (Vector3 &v: my_env->orientation_trail){
            double height = forwardSolver->height_function(v);
            height_log.push_back(height);
        }

        // ambient mesh
        Vector<Vector3> final_rest_positions = my_env->get_new_positions(forwardSolver->inputGeometry->inputVertexPositions.toVector());
        polyscope::registerSurfaceMesh("Bullet resting mesh", final_rest_positions, mesh->getFaceVertexList());
    }
}


int main(int argc, char* argv[])
{
    #ifdef BT_USE_DOUBLE_PRECISION
        printf("BT_USE_DOUBLE_PRECISION\n");
    #else
        printf("Single precision\n");
    #endif

    args::ArgumentParser parser(    "This is a test program.", "This goes after the options.");
    args::HelpFlag help(parser,     "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<std::string> mesh_path_arg(parser, "mesh_path", "path to esh", {'m', "mesh_dir"});
    args::ValueFlag<std::string> center_of_mass_path_arg(parser, "center_of_mass", "center of mass of the shape", {"com"});
	// just get the orientation from 3 input values
	args::ValueFlag<double> o_x_arg(parser, "ox", "orientation x value", {"ox"});
	args::ValueFlag<double> o_y_arg(parser, "oy", "orientation y value", {"oy"});
	args::ValueFlag<double> o_z_arg(parser, "oz", "orientation z value", {"oz"});

    args::Flag drop_flag(parser, "drop", "Perform a quasi-static drop from the given orientation", {"drop"});
    args::Flag probs_flag(parser, "probs", "Compute stable face probabilities and output to file", {"probs"});
    args::ValueFlag<std::string> out_file_arg(parser, "output", "Output file for transformation matrices", {"out"});
    args::Flag viz_flag(parser, "viz", "Visualize orientation and quasi-static trajectory in drop mode", {"viz"});

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
    std::string mesh_path, com_path, out_file;
    if (mesh_path_arg)
        mesh_path = args::get(mesh_path_arg);
    if (center_of_mass_path_arg)
        com_path = args::get(center_of_mass_path_arg);
    if (out_file_arg)
        out_file = args::get(out_file_arg);
    if (o_x_arg && o_y_arg && o_z_arg){
        orientation_gui[0] = args::get(o_x_arg);
        orientation_gui[1] = args::get(o_y_arg);
        orientation_gui[2] = args::get(o_z_arg);
    }

    // extract mesh title
    mesh_title = mesh_path.substr(mesh_path.find_last_of("/\\") + 1);
    mesh_title = mesh_title.substr(0, mesh_title.find_last_of("."));
    // --- Debugging --- //
    polyscope::init();
    vis_utils = VisualUtils();
    load_mesh(mesh_path);
    // Initialize
    update_solver();
    
    // --- Drop mode ---
    if (drop_flag) {
        if (!(o_x_arg && o_y_arg && o_z_arg)) {
            std::cerr << "Error: --ox --oy --oz required for drop mode.\n";
            return 1;
        }
        Vector3 orientation({orientation_gui[0], orientation_gui[1], orientation_gui[2]});
        orientation = orientation.normalize();
        auto [transformation_matrices, trail] = quasi_static_snail_trail_to_global_transformations(forwardSolver, orientation, false);

        if (out_file.empty()) {
            std::cerr << "Error: --out <output_file> required for drop mode.\n";
            return 1;
        }
        // Ensure parent directory exists
        fs::path out_path(out_file);
        if (out_path.has_parent_path()) {
            fs::create_directories(out_path.parent_path());
        }
        std::ofstream ofs(out_file);
        if (!ofs) {
            std::cerr << "Error: Could not open output file: " << out_file << "\n";
            return 1;
        }
        for (const auto& mat : transformation_matrices) {
            ofs << mat << "\n";
        }
        ofs.close();
        std::cout << "Drop transformation matrices written to " << out_file << "\n";
        // return 0;
    }

    // --- Probability mode ---
    if (probs_flag) {
        if (out_file.empty()) {
            std::cerr << "Error: --out <output_file> required for probability mode.\n";
            return 1;
        }
        // Ensure parent directory exists
        fs::path out_path(out_file);
        if (out_path.has_parent_path()) {
            fs::create_directories(out_path.parent_path());
        }
        std::ofstream ofs(out_file);
        if (!ofs) {
            std::cerr << "Error: Could not open output file: " << out_file << "\n";
            return 1;
        }
        ofs << "# face_index normal_x normal_y normal_z probability\n";
        for (Face f : forwardSolver->hullMesh->faces()) {
            if (boundary_builder->face_region_area[f] > 0) {
                Vector3 n = forwardSolver->hullGeometry->faceNormal(f);
                double prob = boundary_builder->face_region_area[f] / (4.0 * M_PI);
                ofs << "Hull face ind: " << f.getIndex() << "  -- normal: "
                    << n.x << " " << n.y << " " << n.z << " \n\t probability: "
                    << prob << "\n";
            }
        }
        ofs.close();
        std::cout << "Stable face probabilities written to " << out_file << "\n";
        // return 0;
    }



    init_visuals(mesh, geometry, forwardSolver, boundary_builder);
    polyscope_defaults();
    polyscope::state::userCallback = myCallback;
    polyscope::show();
    return EXIT_SUCCESS;  
  }