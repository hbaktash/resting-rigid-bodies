// #define BT_USE_DOUBLE_PRECISION

#include "unsupported/Eigen/EulerAngles"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"


#include "nlohmann/json.hpp"
#include "args.hxx"

#include "boundary_tools.h"
#include "mesh_factory.h"
#include "ambient_conversions.h"


// system stuff
#include "chrono"
#include "file_IO.h"

// Polyscope (optional)
#ifdef RESTING_RIGID_BODIES_USE_POLYSCOPE
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "imgui.h"
#include "visual_utils.h"
#endif

// Bullet (optional)
#ifdef RESTING_RIGID_BODIES_USE_BULLET
#include "LinearMath/btVector3.h"
#include "btBulletDynamicsCommon.h"
#include "bullet_sim.h"
#endif


using namespace geometrycentral;
using namespace geometrycentral::surface;


// == Global vars
// QS simulation and probabilities
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;
Vector3 G;
Forward3DSolver* forwardSolver;
BoundaryBuilder *boundary_builder;
float orientation_gui[3] = {0., 0., 1.};
double min_QS_rotation_angle = 1.; // 0.005 for slide videos. 1 means no refinement (contact events only)
Vector3 down_vec = Vector3({0, -1, 0});

std::string out_file;

// DEBUG visuals and Polyscope globals (guarded)
#ifdef RESTING_RIGID_BODIES_USE_POLYSCOPE
int drop_step_slider = 0;
std::vector<VertexData<Vector3>> drop_step_positions;
std::vector<Vector3> drop_step_orientations;
VisualUtils vis_utils;

static void polyscope_defaults() {
  polyscope::options::ssaaFactor = 2;
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  if (polyscope::hasCurveNetwork("Arc curves all edge arcs"))
    polyscope::getCurveNetwork("Arc curves all edge arcs")->setEnabled(false);
  if (polyscope::hasPointCloud("Edge equilibria"))
    polyscope::getPointCloud("Edge equilibria")->setEnabled(false);
  if (polyscope::hasPointCloud("stable Vertices Normals"))
    polyscope::getPointCloud("stable Vertices Normals")->setEnabled(false);
  if (polyscope::hasPointCloud("Center of Mass"))
    polyscope::getPointCloud("Center of Mass")->setEnabled(false);
}
#endif

// Bullet globals and helpers (guarded)
#ifdef RESTING_RIGID_BODIES_USE_BULLET
PhysicsEnv* my_env = nullptr;
bool save_Bullet_trail = false;

static void initialize_bullet_env(Vector3 orientation){
  double bullet_step_size = 0.01; // 0.01 for slide videos
  double ground_box_y = -2.1;
  Vector3 ground_box_shape({10,1,10});
  my_env = new PhysicsEnv();
  my_env->init_physics();
  my_env->init_geometry(forwardSolver->hullMesh, forwardSolver->hullGeometry);
  my_env->add_ground(ground_box_y, ground_box_shape);
  my_env->add_object(G, orientation);
  my_env->default_step_size = bullet_step_size;
}
#endif

// Polyscope callback (guarded)
#ifdef RESTING_RIGID_BODIES_USE_POLYSCOPE
void myCallback() {
  // sliders for orientation
  if (ImGui::SliderFloat3("orientation", orientation_gui, -10.0f, 10.0f)) {
    Vector3 orientation = Vector3({orientation_gui[0], orientation_gui[1], orientation_gui[2]}).normalize();
    std::vector<Vector3> orientation_vec = {orientation * vis_utils.gm_radi + vis_utils.center};
    polyscope::registerPointCloud("orientation on gm", orientation_vec)
      ->setPointRadius(0.03, false)->setPointColor({0,0,0});
  }
  if (ImGui::Button("show orientation and QS trail")) {
    Vector3 orientation({orientation_gui[0], orientation_gui[1], orientation_gui[2]});
    orientation = orientation.normalize();
    Vector3 down_vec({0,-1,0});
    VertexData<Vector3> new_positions = orient_mesh_to_unit_vec(orientation, down_vec, mesh, geometry);
    polyscope::registerSurfaceMesh("oriented mesh", new_positions, mesh->getFaceVertexList());
    std::vector<Vector3> orientation_vec = {orientation * vis_utils.gm_radi + vis_utils.center};
    polyscope::registerPointCloud("orientation on gm", orientation_vec)->setPointRadius(0.03, false)->setPointColor({0,0,0});
    std::vector<Vector3> snail_trail = forwardSolver->snail_trail_log(orientation);
    draw_trail_on_gm(snail_trail, {39./255., 189./255., 0}, "quasi-static trail", 3.);
  }
  if (ImGui::Button("save orientation trajectory to file")) {
    if (out_file.empty()) {
      std::cerr << "Error: --out <output_file> required for saving.\n";
    } else {
      Vector3 orientation({orientation_gui[0], orientation_gui[1], orientation_gui[2]});
      orientation = orientation.normalize();
      auto [transformation_matrices, saved_snail_trail_refined] =
          QS_trail_to_global_transformations(forwardSolver, orientation, down_vec, min_QS_rotation_angle);
      save_trans_mats_and_orientations_to_file(transformation_matrices, saved_snail_trail_refined, out_file);
    }
  }
  if (ImGui::Button("show drop sequence")) {
    //hide initial input mesh and hull, and the oriented mesh
    polyscope::getSurfaceMesh("oriented mesh")->setEnabled(false);
    polyscope::getSurfaceMesh("init hull mesh")->setEnabled(false);
    polyscope::getSurfaceMesh("init input mesh")->setEnabled(false);
    // save drop positions and show the first one
    Vector3 orientation({orientation_gui[0], orientation_gui[1], orientation_gui[2]});
    orientation = orientation.normalize();
    auto [transformation_matrices, saved_snail_trail_refined] =
        QS_trail_to_global_transformations(forwardSolver, orientation, down_vec, min_QS_rotation_angle);
    save_trans_mats_and_orientations_to_file(transformation_matrices, saved_snail_trail_refined, out_file);
    draw_ground_plane_mesh(down_vec, 0);
    drop_step_positions.clear();
    drop_step_orientations.clear();
    for (int i = 0; i < (int)transformation_matrices.size(); i++) {
      VertexData<Vector3> new_positions = apply_trans_to_positions(geometry->inputVertexPositions, transformation_matrices[i]);
      drop_step_positions.push_back(new_positions);
      drop_step_orientations.push_back(saved_snail_trail_refined[i]);
    }
    // show the first step
    Vector3 orientation0 = drop_step_orientations[0];
    polyscope::registerSurfaceMesh("drop step", drop_step_positions[0], mesh->getFaceVertexList());
    std::vector<Vector3> orientation_vec = {orientation0 * vis_utils.gm_radi + vis_utils.center};
    polyscope::registerPointCloud("orientation on gm", orientation_vec)->setPointRadius(0.03, false)->setPointColor({0,0,0});
  }
  if (ImGui::SliderInt("drop step", &drop_step_slider, 0, (int)drop_step_orientations.size()-1)) {
    Vector3 orientation = drop_step_orientations[drop_step_slider];
    polyscope::registerSurfaceMesh("drop step", drop_step_positions[drop_step_slider], mesh->getFaceVertexList());
    std::vector<Vector3> orientation_vec = {orientation * vis_utils.gm_radi + vis_utils.center};
    polyscope::registerPointCloud("orientation on gm", orientation_vec)->setPointRadius(0.03, false)->setPointColor({0,0,0});
  }

  #ifdef RESTING_RIGID_BODIES_USE_BULLET
  ImGui::Checkbox("save bullet orientation trail", &save_Bullet_trail);
  if (ImGui::Button("bullet snail trail")) {
    Vector3 orientation({orientation_gui[0], orientation_gui[1], orientation_gui[2]});
    orientation = orientation.normalize();
    initialize_bullet_env(orientation, false);
    Face touching_face = my_env->final_stable_face(true);
    std::vector<Vector3> snail_trail = my_env->orientation_trail;
    for (Vector3 &v: snail_trail) v += vis_utils.center;
    auto trail_pc = polyscope::registerPointCloud("bullet pc trail", snail_trail);
    glm::vec3 bullet_color = {209./255., 227./255., 28./255.};
    trail_pc->setPointColor(bullet_color);
    trail_pc->setPointRadius(0.003, false);
    polyscope::registerPointCloud("final bullet orientation", std::vector<Vector3>{snail_trail.back()})
      ->setPointColor({1,0,0})->setPointRadius(0.03, false);
    std::vector<double> height_log;
    for (Vector3 &v: my_env->orientation_trail) height_log.push_back(forwardSolver->height_function(v));
    Vector<Vector3> final_rest_positions = my_env->get_new_positions(forwardSolver->inputGeometry->inputVertexPositions.toVector());
    polyscope::registerSurfaceMesh("Bullet resting mesh", final_rest_positions, mesh->getFaceVertexList());
  }
  #endif
}
#endif // RESTING_RIGID_BODIES_USE_POLYSCOPE

void load_mesh(std::string path, bool preprocess = true){
	std::unique_ptr<SurfaceMesh> nm_mesh_ptr;
    std::unique_ptr<VertexPositionGeometry> geometry_ptr;
	std::unique_ptr<VertexPositionGeometry> nm_geometry_ptr;
	std::tie(nm_mesh_ptr, nm_geometry_ptr) = readSurfaceMesh(path);
	SurfaceMesh *nm_mesh = nm_mesh_ptr.release();
	VertexPositionGeometry *nm_geometry = nm_geometry_ptr.release();
	
	nm_mesh->greedilyOrientFaces();
	nm_mesh->compress();
    std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
	mesh_ptr = nm_mesh->toManifoldMesh();
	mesh = mesh_ptr.release();
	geometry = new VertexPositionGeometry(*mesh);
	// transfer from nm geometry
	for (Vertex v : mesh->vertices()) {
		geometry->inputVertexPositions[v.getIndex()] = nm_geometry->inputVertexPositions[v.getIndex()];
	}

	// preproccess and shift COM
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


int main(int argc, char* argv[])
{
    #ifdef BT_USE_DOUBLE_PRECISION
        printf("BT_USE_DOUBLE_PRECISION\n");
    #else
        printf("Single precision\n");
    #endif

    args::ArgumentParser parser(    "This is a test program.", "This goes after the options.");
    args::HelpFlag help(parser,     "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<std::string> mesh_path_arg(parser, "mesh_path", "path to esh", {'m', "mesh"});
    args::ValueFlag<std::string> center_of_mass_path_arg(parser, "center_of_mass", "center of mass of the shape", {"com"});
    // just get the orientation from 3 input values
    args::ValueFlag<double> o_x_arg(parser, "ox", "orientation x value", {"ox"});
    args::ValueFlag<double> o_y_arg(parser, "oy", "orientation y value", {"oy"});
    args::ValueFlag<double> o_z_arg(parser, "oz", "orientation z value", {"oz"});

    args::Flag drop_flag(parser, "drop", "Perform a quasi-static drop from the given orientation", {"drop"});
    args::ValueFlag<double> min_QS_rotation_angle_arg(parser, "min_QS_rotation_angle", "Minimum quasi-static rotation angle", {"min_QS_angle"});
    args::Flag probs_flag(parser, "probs", "Compute stable face probabilities and output to file", {"probs"});
    args::ValueFlag<std::string> out_file_arg(parser, "output", "Output file for transformation matrices / probabilities", {"out"});
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
    std::string mesh_path, com_path;
    if (mesh_path_arg)
        mesh_path = args::get(mesh_path_arg);
    if (center_of_mass_path_arg)
        com_path = args::get(center_of_mass_path_arg);
    if (out_file_arg)
        out_file = args::get(out_file_arg);
    if (min_QS_rotation_angle_arg)
        min_QS_rotation_angle = args::get(min_QS_rotation_angle_arg);
    if (o_x_arg && o_y_arg && o_z_arg){
        orientation_gui[0] = args::get(o_x_arg);
        orientation_gui[1] = args::get(o_y_arg);
        orientation_gui[2] = args::get(o_z_arg);
    }

    // Load mesh and initialize solver (headless-friendly)
    load_mesh(mesh_path);
    update_solver();

    // Require output file for headless outputs
    if ((drop_flag || probs_flag) && out_file.empty() && !viz_flag) {
        std::cerr << "Error: --out <output_file> required for saving. Or use visualize mode\n";
        return 1;
    }

    // --- Drop mode ---
    if (drop_flag) {
        if (!(o_x_arg && o_y_arg && o_z_arg)) {
            std::cerr << "Error: --ox --oy --oz required for drop mode.\n";
            return 1;
        }
        Vector3 orientation({orientation_gui[0], orientation_gui[1], orientation_gui[2]});
        orientation = orientation.normalize();

        auto [transformation_matrices, refined_trail] =
            QS_trail_to_global_transformations(forwardSolver, orientation, down_vec, min_QS_rotation_angle);

        // Ensure parent dir exists and write JSON/txt
        save_trans_mats_and_orientations_to_file(transformation_matrices, refined_trail, out_file);

        // Optional visualization in drop mode
        if (viz_flag) {
        #ifdef RESTING_RIGID_BODIES_USE_POLYSCOPE
            polyscope::init();
            init_visuals(mesh, geometry, forwardSolver, boundary_builder);
            polyscope_defaults();

            // Input orientation on GM
            std::vector<Vector3> orientation_vec = {orientation * vis_utils.gm_radi + vis_utils.center};
            polyscope::registerPointCloud("input orientation on gm", orientation_vec)
              ->setPointRadius(0.03, false)->setPointColor({0,0,0});

            // Quasi-static trail on GM
            std::vector<Vector3> snail_trail = forwardSolver->snail_trail_log(orientation);
            draw_trail_on_gm(snail_trail, {39./255., 189./255., 0}, "quasi-static trail", 3.);
            polyscope::show();
        #else
            std::cerr << "Warning: --viz requested but Polyscope was disabled at build time.\n";
        #endif
        }
        return EXIT_SUCCESS;
    }

    // --- Probability mode ---
    if (probs_flag) {
        // already computed in update_solver()
        Eigen::VectorXd face_region_areas = boundary_builder->face_region_area.toVector();
        forwardSolver->hullGeometry->requireFaceNormals();
        Eigen::VectorX<Vector3> face_normals = forwardSolver->hullGeometry->faceNormals.toVector();
        forwardSolver->hullGeometry->unrequireFaceNormals();
        std::vector<Eigen::Matrix4d> stable_trans_mats;
        for (Face f: forwardSolver->hullMesh->faces()){
            if (face_region_areas[f.getIndex()] > 0) {
                Eigen::Matrix4d trans_mat = trans_mat_for_orientation(
                    face_normals[f.getIndex()], down_vec, forwardSolver->hullMesh, forwardSolver->hullGeometry);
                stable_trans_mats.push_back(trans_mat);
            }
        }
        save_probabilities_to_file(face_region_areas, face_normals, stable_trans_mats, out_file);
    }

    // --- Visualization mode (interactive) ---
    if (viz_flag) {
    #ifdef RESTING_RIGID_BODIES_USE_POLYSCOPE
        polyscope::init();
        init_visuals(mesh, geometry, forwardSolver, boundary_builder);
        polyscope_defaults();
        polyscope::state::userCallback = myCallback;
        polyscope::show();
    #else
        std::cerr << "Warning: --viz requested but Polyscope was disabled at build time.\n";
    #endif
    }

    return EXIT_SUCCESS;
  }