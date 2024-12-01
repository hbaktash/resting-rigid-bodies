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
#include "deformation.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


#define ANSI_FG_MAGENTA "\x1b[35m"
#define ANSI_FG_YELLOW "\x1b[33m"
#define ANSI_FG_GREEN "\x1b[32m"
#define ANSI_FG_WHITE "\x1b[37m"
#define ANSI_FG_RED "\x1b[31m"
#define ANSI_RESET "\x1b[0m"


// == Geometry-central data
ManifoldSurfaceMesh *hull_mesh, *concave_mesh, *optimized_concave_mesh;
VertexPositionGeometry *hull_geometry, *concave_geometry, *optimized_concave_geometry;
Vector3 goal_G, // center of Mass of pre-deformed hull
        current_G;

bool load_hull_from_file = true;

Forward3DSolver* forwardSolver;
  

VisualUtils vis_utils;

ManifoldSurfaceMesh* sphere_mesh;
VertexPositionGeometry* sphere_geometry;

// raster image stuff
FaceData<Vector3> face_colors;

// boundary stuff
BoundaryBuilder *boundary_builder;

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
      cluster_distance_reg = 0.0,
      unstable_attraction_thresh = 0.1;
int sobolev_p = 2;
// optimization stuff



// example choice
std::string hull_input_name, concave_input_name;
    

void visualize_gauss_map(Forward3DSolver* forwardSolver){
  std::unique_ptr<ManifoldSurfaceMesh> sphere_mesh_ptr;
  std::unique_ptr<VertexPositionGeometry> sphere_geometry_ptr;
  std::tie(sphere_mesh_ptr, sphere_geometry_ptr) = generate_polyhedra("sphere");
  sphere_mesh = sphere_mesh_ptr.release();
  sphere_geometry = sphere_geometry_ptr.release();
  
  // update vis utils
  vis_utils.draw_gauss_map(forwardSolver, sphere_mesh, sphere_geometry);
}


void update_solver(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry, Vector3 G = Vector3({-1,-1,-1})){ // only doing this for convex input
  forwardSolver = new Forward3DSolver(mesh, geometry, G, true);
  if (G == Vector3({-1,-1,-1})){ // not provided
    forwardSolver->set_uniform_G();
  }
  forwardSolver->initialize_pre_computes();
  boundary_builder = new BoundaryBuilder(forwardSolver);
  boundary_builder->build_boundary_normals();
}




void initialize_state(std::string hull_input_name, std::string concave_input_name){
    std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
    std::unique_ptr<VertexPositionGeometry> geometry_ptr;
    
    std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(hull_input_name);
    hull_mesh = mesh_ptr.release();
    hull_geometry = geometry_ptr.release();

    std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(concave_input_name);
    concave_mesh = mesh_ptr.release();
    concave_geometry = geometry_ptr.release();
    auto currentGVpair = find_center_of_mass(*concave_mesh, *concave_geometry);
    current_G = currentGVpair.first;

    bool triangulate = true;
    // preprocess_mesh(hull_mesh, hull_geometry, triangulate, false, 1.);
    preprocess_mesh(concave_mesh, concave_geometry, triangulate, false, 1.);

    // no need for path anymore
    hull_input_name = hull_input_name.substr(hull_input_name.find_last_of("/") + 1);
    hull_input_name = hull_input_name.substr(0, hull_input_name.find_last_of("."));
    concave_input_name = concave_input_name.substr(concave_input_name.find_last_of("/") + 1);
    concave_input_name = concave_input_name.substr(0, concave_input_name.find_last_of("."));


    update_solver(hull_mesh, hull_geometry, goal_G);

    // visuals
    polyscope::registerSurfaceMesh("hull input",
                                hull_geometry->inputVertexPositions, hull_mesh->getFaceVertexList())->setTransparency(0.7);
    polyscope::registerSurfaceMesh("concave input",
                                concave_geometry->inputVertexPositions, concave_mesh->getFaceVertexList())->setTransparency(0.7); 
    polyscope::registerPointCloud("goal G", std::vector<Vector3>{goal_G})->setPointColor({0,0,0})->setPointRadius(vis_utils.G_radi);
    polyscope::registerPointCloud("current G", std::vector<Vector3>{current_G})->setPointColor({1,0,0})->setPointRadius(vis_utils.G_radi);
}


// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
  ImGui::Checkbox("adaptive reg", &adaptive_reg);
  ImGui::Checkbox("visualize steps", &visualize_steps);
  if (ImGui::Button("save optimized shape & G")){
    std::string output_name = "filled_in_" + concave_input_name;
    writeSurfaceMesh(*optimized_concave_mesh, *optimized_concave_geometry, "../meshes/hulls/" + std::string(output_name) +".obj");
    // save G using filesystem
    std::ofstream out_file("../meshes/hulls/" + std::string(output_name) + "_finalG.txt");
    out_file << goal_G.x << " " << goal_G.y << " " << goal_G.z << "\n";
    out_file.close();
  }
}


int main(int argc, char **argv) {
    // Parse args
    args::ArgumentParser parser("Inverse convex hull filling");
    // make parser for input name and other params
    args::HelpFlag help(parser,    "help", "Display this help menu", {'h', "help"});
    // args::Flag do_just_ours(parser, "do_just_ours", "do just ours", {"just_ours"});
    // args::ValueFlag<int> total_samples(parser, "ICOS_samples", "Total number of samples", {"samples"});
    args::ValueFlag<std::string> convex_input_arg(parser, "convex_input_str", "path to input convex shape to be filled", {'c', "convex_dir"});
    args::ValueFlag<std::string> concave_input_arg(parser, "concave_input_arg", "path to concave shape to fill with", {'m', "concave_dir"});
    args::ValueFlag<std::string> goal_G_text_arg(parser, "goal_G_text_arg", "path to text containing G", {'g', "com"});

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

    hull_input_name  = args::get(convex_input_arg);
    concave_input_name = args::get(concave_input_arg);
    std::string goal_G_text_name = args::get(goal_G_text_arg);
    std::ifstream in_file(goal_G_text_name);
    in_file >> goal_G.x >> goal_G.y >> goal_G.z;
    in_file.close();

    std::cout << ANSI_FG_YELLOW << "input convex shape: " << hull_input_name << ANSI_RESET << "\n";
    std::cout << ANSI_FG_YELLOW << "input concave shape: " << concave_input_name << ANSI_RESET << "\n";
    std::cout << ANSI_FG_YELLOW << "goal G: " << goal_G << ANSI_RESET << "\n";

    // build mesh
    vis_utils = VisualUtils();
    polyscope::init();

    initialize_state(hull_input_name, concave_input_name);
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