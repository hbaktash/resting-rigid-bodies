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
ManifoldSurfaceMesh *hull_mesh, *concave_mesh, *optimized_concave_mesh, *ref_remeshed_concave_mesh;
VertexPositionGeometry *hull_geometry, *concave_geometry, *optimized_concave_geometry, *ref_remeshed_concave_geometry;
Vector3 goal_G, // center of Mass of pre-deformed hull
        current_G;

VisualUtils vis_utils;

ManifoldSurfaceMesh* sphere_mesh;
VertexPositionGeometry* sphere_geometry;

// raster image stuff
FaceData<Vector3> face_colors;



// deformation stuff
bool animate_deform = false,
     animate_G_deform = false,
     v2_dice_animate = false,
     enforce_snapping = false,
     use_reg = false,
     use_QP_solver = true,
     curvature_weighted_CP = false,
     dynamic_remesh = true;

float dice_search_decay = 0.95;

float bending_lambda_exps[2] = {1., 1.},
      membrane_lambda_exps[2] = {3., 3.},
      CP_lambda_exps[2] = {1., 7.},
      barrier_lambda_exps[2] = {-4., -8.},
      G_lambda_exps[2] = {1,5},
      reg_lambda_exp = -3.,
      internal_p = 0.91,
      refinement_CP_threshold = 0.001,
      active_set_threshold = 0.01,
      split_robustness_threshold = 0.2;
int filling_max_iter = 10;

// deformed shape per step
std::vector<std::pair<ManifoldSurfaceMesh*, VertexPositionGeometry*>> deformed_shapes;


// example choice
std::string hull_input_name, concave_input_name,
            hull_shape_name, concave_shape_name;
    

void visualize_gauss_map(Forward3DSolver* forwardSolver){
  std::unique_ptr<ManifoldSurfaceMesh> sphere_mesh_ptr;
  std::unique_ptr<VertexPositionGeometry> sphere_geometry_ptr;
  std::tie(sphere_mesh_ptr, sphere_geometry_ptr) = generate_polyhedra("sphere");
  sphere_mesh = sphere_mesh_ptr.release();
  sphere_geometry = sphere_geometry_ptr.release();
  
  // update vis utils
  vis_utils.draw_gauss_map(forwardSolver, sphere_mesh, sphere_geometry);
}


void ideal_hull_G_stats(bool visualize = true){
    // ideal G and hull probs
    Forward3DSolver* ideal_solver = new Forward3DSolver(hull_mesh, hull_geometry, goal_G, false); //
    ideal_solver->initialize_pre_computes();
    BoundaryBuilder* ideal_bnd_builder = new BoundaryBuilder(ideal_solver);
    ideal_bnd_builder->build_boundary_normals();
    std::cout << "ideal G:" << goal_G << "\n";
    std::cout << "ideal probabilities:\n";
    ideal_bnd_builder->print_area_of_boundary_loops();
    if (visualize){
        visualize_gauss_map(ideal_solver);
        vis_utils.update_visuals(ideal_solver, ideal_bnd_builder, sphere_mesh, sphere_geometry);
        vis_utils.draw_stable_patches_on_gauss_map(false, ideal_bnd_builder, false);
    }
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
    // preprocess_mesh(concave_mesh, concave_geometry, triangulate, false, 1.);

    // no need for path anymore
    hull_shape_name = hull_input_name.substr(hull_input_name.find_last_of("/") + 1);
    hull_shape_name = hull_shape_name.substr(0, hull_shape_name.find_last_of("."));
    concave_shape_name = concave_input_name.substr(concave_input_name.find_last_of("/") + 1);
    concave_shape_name = concave_shape_name.substr(0, concave_shape_name.find_last_of("."));

    // ideal stats
    ideal_hull_G_stats(false);

    // visuals
    polyscope::registerSurfaceMesh("hull input",
                                hull_geometry->inputVertexPositions, hull_mesh->getFaceVertexList())->setTransparency(0.7);
    polyscope::registerSurfaceMesh("concave input",
                                concave_geometry->inputVertexPositions, concave_mesh->getFaceVertexList())->setTransparency(0.7); 
    polyscope::registerPointCloud("goal G", std::vector<Vector3>{goal_G})->setPointColor({0,0,0})->setPointRadius(vis_utils.G_radi);
    polyscope::registerPointCloud("current G", std::vector<Vector3>{current_G})->setPointColor({1,0,0})->setPointRadius(vis_utils.G_radi);
}


void initialize_deformation_params(DeformationSolver *deformation_solver){
  deformation_solver->dynamic_remesh = dynamic_remesh;
  deformation_solver->filling_max_iter = filling_max_iter;  
  
  deformation_solver->init_bending_lambda = pow(10, bending_lambda_exps[0]);
  deformation_solver->final_bending_lambda = pow(10, bending_lambda_exps[1]);
  deformation_solver->init_membrane_lambda = pow(10, membrane_lambda_exps[0]);
  deformation_solver->final_membrane_lambda = pow(10, membrane_lambda_exps[1]);
  deformation_solver->init_CP_lambda = pow(10, CP_lambda_exps[0]);
  deformation_solver->final_CP_lambda = pow(10, CP_lambda_exps[1]);
  if (!use_QP_solver){
    deformation_solver->init_barrier_lambda = pow(10, barrier_lambda_exps[0]);
    deformation_solver->final_barrier_lambda = pow(10, barrier_lambda_exps[1]);
  }
  else {
    deformation_solver->init_barrier_lambda = 0.;
    deformation_solver->final_barrier_lambda = 0.;
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


// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {

    // deformation params
    ImGui::Checkbox("do remeshing", &dynamic_remesh);
    ImGui::Checkbox("enforce snapping at threshold", &enforce_snapping);
    // ImGui::Checkbox("curvature weighted", &curvature_weighted_CP);
    ImGui::SliderInt("filling iters", &filling_max_iter, 1, 300);
    ImGui::SliderFloat("growth p", &internal_p, 0., 1.);
    ImGui::SliderFloat("refinement CP threshold ", &refinement_CP_threshold, 0., 1.);
    ImGui::SliderFloat("split rubostness threshold ", &split_robustness_threshold, 0., 1.);
    ImGui::InputFloat2("init/final bending log ", bending_lambda_exps);
    ImGui::InputFloat2("init/final membrane log ", membrane_lambda_exps);
    ImGui::InputFloat2("init/final CP log ", CP_lambda_exps);
    ImGui::Checkbox("use QP solver ", &use_QP_solver);
    if(!use_QP_solver) 
        ImGui::InputFloat2("init/final barrier log ", barrier_lambda_exps);
    else
        ImGui::InputFloat("active set threshold ", &active_set_threshold);
    ImGui::Checkbox("use reg ", &use_reg);
    if(use_reg) 
        ImGui::SliderFloat("reg lambda; log10", &reg_lambda_exp, -5., 5.);
    

    if (ImGui::Button("deform to hull")){
        animate_deform = true;
    }

    if (ImGui::Button("refresh")){
        polyscope::removeAllStructures();
        initialize_state(hull_input_name, concave_input_name);
    }

    if (ImGui::Button("save optimized concave shape")){
        std::string output_name = "hull_fill_output_" + concave_shape_name;
        writeSurfaceMesh(*optimized_concave_mesh, *optimized_concave_geometry, "../meshes/hulls/" + std::string(output_name) +".obj");
        std::string ref_output_name = "hull_fill_remeshed_ref_" + concave_shape_name;
        writeSurfaceMesh(*ref_remeshed_concave_mesh, *ref_remeshed_concave_geometry, "../meshes/hulls/" + std::string(ref_output_name) +".obj");
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
    while (true) {
        if (animate_deform){
            animate_deform = false;
            DeformationSolver* deformationSolver = new DeformationSolver(concave_mesh, concave_geometry, hull_mesh, hull_geometry);   
            initialize_deformation_params(deformationSolver);
            size_t current_fill_iter = 0;
            Eigen::MatrixXd new_points = deformationSolver->solve_for_bending(1, false, nullptr, nullptr);

            optimized_concave_mesh = deformationSolver->mesh;
            optimized_concave_geometry = deformationSolver->deformed_geometry;
            ref_remeshed_concave_mesh = deformationSolver->mesh;
            ref_remeshed_concave_geometry = deformationSolver->old_geometry;

            // checking probabilities for the deformed shape
            Forward3DSolver* final_solver = new Forward3DSolver(deformationSolver->mesh, deformationSolver->deformed_geometry, goal_G, true); // 
            
            polyscope::registerSurfaceMesh("deformed final hull", final_solver->hullGeometry->inputVertexPositions, 
                                           final_solver->hullMesh->getFaceVertexList())->setTransparency(0.4);

            
            // post process
            // get the probabilities of new deformed shape
            final_solver->set_uniform_G();
            final_solver->initialize_pre_computes();
            BoundaryBuilder* tmp_bnd_builder = new BoundaryBuilder(final_solver);
            tmp_bnd_builder->build_boundary_normals();

            Vector3 post_deform_G = final_solver->get_G();
            std::cout << "post-deform G:" << post_deform_G << "\n";
            printf(" post deform probabilities:\n");
            tmp_bnd_builder->print_area_of_boundary_loops();
            // visuals for post deform
            // None yet

            // ideal stats
            ideal_hull_G_stats(false);
        }
        polyscope::frameTick();
    }
    polyscope::show();

    return EXIT_SUCCESS;
}


//