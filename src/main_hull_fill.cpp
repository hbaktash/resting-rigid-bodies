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
ManifoldSurfaceMesh *hull_mesh, *concave_mesh, 
                    *optimized_concave_mesh, *ref_remeshed_concave_mesh,
                    *deformed_mesh;
VertexPositionGeometry *hull_geometry, *concave_geometry, 
                       *optimized_concave_geometry, *ref_remeshed_concave_geometry,
                       *deformed_geometry;
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
     dynamic_remesh = false,
     use_static_dists = true;

double concave_input_remesh_scale = 1.;

float bending_lambda_exps[2]  = {0. , 0.},
      membrane_lambda_exps[2] = {0.5, 0.5},
      CP_lambda_exps[2]       = {2. , 7.},
      barrier_lambda_exps[2]  = {-4., -8.},
      G_lambda_exps[2]        = {2  , 7},
      reg_lambda_exp          = 1.,
      G_deform_sobolev_lambda = 10.,
      internal_p              = 0.2,
      refinement_CP_threshold = 0.00,
      active_set_threshold    = 0.01,
      split_robustness_threshold = 0.0;
int filling_max_iter = 100;

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
    std::cout << "ideal G stats G\n \t\t:" << goal_G << "\n";
    std::cout << "ideal probabilities:\n";
    ideal_bnd_builder->print_area_of_boundary_loops();
    if (visualize){
        visualize_gauss_map(ideal_solver);
        vis_utils.update_visuals(ideal_solver, ideal_bnd_builder, sphere_mesh, sphere_geometry);
        vis_utils.draw_stable_patches_on_gauss_map(false, ideal_bnd_builder, false);
        
        Vector3 offset({2., 0., 0});
        auto ps_shifted_hull = polyscope::registerSurfaceMesh("offset hull", ideal_solver->hullGeometry->inputVertexPositions + offset, ideal_solver->hullMesh->getFaceVertexList());
        
        FaceData<double> probs = ideal_bnd_builder->face_region_area/(4.*PI);
        ps_shifted_hull->addFaceScalarQuantity("raw probs ", probs)->setColorMap("reds")->setEnabled(true);
        for (Face f: hull_mesh->faces()){
            probs[f] = probs[ideal_solver->face_last_face[f]];
        }
        ps_shifted_hull->addFaceScalarQuantity("probs accum", probs)->setColorMap("reds")->setEnabled(true);
    }
}


void initialize_state(std::string hull_input_name, std::string concave_input_name, bool normalize_concave_input = false){
    std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
    std::unique_ptr<VertexPositionGeometry> geometry_ptr;
    
    std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(hull_input_name);
    hull_mesh = mesh_ptr.release();
    hull_geometry = geometry_ptr.release();

    std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(concave_input_name);
    concave_mesh = mesh_ptr.release();
    concave_geometry = geometry_ptr.release();
    bool triangulate = true;
    // preprocess_mesh(hull_mesh, hull_geometry, triangulate, false, 1.);
    if (normalize_concave_input){
        preprocess_mesh(concave_mesh, concave_geometry, triangulate, concave_input_remesh_scale != 1. ? true : false, concave_input_remesh_scale);
    }

    auto currentGVpair = find_center_of_mass(*concave_mesh, *concave_geometry);
    current_G = currentGVpair.first;


    // no need for path anymore
    hull_shape_name = hull_input_name.substr(hull_input_name.find_last_of("/") + 1);
    hull_shape_name = hull_shape_name.substr(0, hull_shape_name.find_last_of("."));
    concave_shape_name = concave_input_name.substr(concave_input_name.find_last_of("/") + 1);
    concave_shape_name = concave_shape_name.substr(0, concave_shape_name.find_last_of("."));


    // visuals
    polyscope::registerSurfaceMesh("hull input",
                                hull_geometry->inputVertexPositions, hull_mesh->getFaceVertexList())->setTransparency(0.7);
    polyscope::registerSurfaceMesh("concave input",
                                concave_geometry->inputVertexPositions, concave_mesh->getFaceVertexList())->setTransparency(0.7); 
    polyscope::registerPointCloud("goal G", std::vector<Vector3>{goal_G})->setPointColor({0,0,0})->setPointRadius(0.01);
    polyscope::registerPointCloud("current G", std::vector<Vector3>{current_G})->setPointColor({1,0,0})->setPointRadius(0.01);
    
    // ideal stats
    ideal_hull_G_stats(false);
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
    deformation_solver->init_G_lambda = pow(10, G_lambda_exps[0]);
    deformation_solver->final_G_lambda = pow(10, G_lambda_exps[1]);
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
    deformation_solver->G_deform_sobolev_lambda = G_deform_sobolev_lambda;
    deformation_solver->use_static_dists = use_static_dists;
    deformation_solver->use_QP_solver = use_QP_solver;
}

void save_params_to_file(std::string file_path){
    std::ofstream out_file(file_path);
    out_file << "dynamic_remesh: " << dynamic_remesh << std::endl;
    out_file << "filling_max_iter " << filling_max_iter << "\n";
    out_file << "enforce_snapping: " << enforce_snapping << std::endl;
    out_file << "curvature_weighted_CP: " << curvature_weighted_CP << std::endl;
    out_file << "filling_max_iter: " << filling_max_iter << std::endl;
    out_file << "internal_p: " << internal_p << std::endl;
    out_file << "refinement_CP_threshold: " << refinement_CP_threshold << std::endl;
    out_file << "split_robustness_threshold: " << split_robustness_threshold << std::endl;
    out_file << "bending_lambda_exps: " << bending_lambda_exps[0] << " " << bending_lambda_exps[1] << std::endl;
    out_file << "membrane_lambda_exps: " << membrane_lambda_exps[0] << " " << membrane_lambda_exps[1] << std::endl;
    out_file << "CP_lambda_exps: " << CP_lambda_exps[0] << " " << CP_lambda_exps[1] << std::endl;
    out_file << "G_lambda_exps: " << G_lambda_exps[0] << " " << G_lambda_exps[1] << std::endl;
    out_file << "use_QP_solver: " << use_QP_solver << std::endl;
    out_file << "active_set_threshold: " << active_set_threshold << std::endl;
    out_file << "use_reg: " << use_reg << std::endl;
    out_file << "reg_lambda_exp: " << reg_lambda_exp << std::endl;
    out_file << "use_static_dists for G-def" << use_static_dists << "\n";
    out_file.close();
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {

    // deformation params
    ImGui::Checkbox("dynamic remeshing", &dynamic_remesh);
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

    ImGui::InputFloat2("init/final G log ", G_lambda_exps);
    ImGui::InputFloat("G defrom sobolev lambda ", &G_deform_sobolev_lambda);
    ImGui::Checkbox("use static dists for G-deform ", &use_static_dists);
    if (ImGui::Button("do the G-deform")){
        animate_G_deform = true;
    }

    if (ImGui::Button("refresh")){
        polyscope::removeAllStructures();
        initialize_state(hull_input_name, concave_input_name);
    }

    if (ImGui::Button("save optimized concave shape")){
        std::string output_name = "hull_fill_output_" + concave_shape_name;
        writeSurfaceMesh(*optimized_concave_mesh, *optimized_concave_geometry, "../meshes/hulls/nonconvex_deformation/deformation_saves/" 
                                                                            + std::string(output_name) +".obj");
        std::string ref_output_name = "hull_fill_remeshed_ref_" + concave_shape_name;
        writeSurfaceMesh(*ref_remeshed_concave_mesh, *ref_remeshed_concave_geometry, "../meshes/hulls/nonconvex_deformation/deformation_saves/" 
                                                                            + std::string(ref_output_name) +".obj");
        // save parameters in a text file
        std::string output_params_path = "../meshes/hulls/nonconvex_deformation/deformation_saves/params_" + std::string(output_name) + ".txt";
        save_params_to_file(output_params_path);
    }
    if (ImGui::Button("save optimized G-deform shape")){
        std::string output_name = "G_deform_output_" + concave_shape_name;
        writeSurfaceMesh(*optimized_concave_mesh, *optimized_concave_geometry, "../meshes/hulls/nonconvex_deformation/deformation_saves/" 
                                                                            + std::string(output_name) +".obj");
        // save parameters in a text file
        std::string output_params_path = "../meshes/hulls/nonconvex_deformation/deformation_saves/params_" + std::string(output_name) + ".txt";
        save_params_to_file(output_params_path);
    }
}


int main(int argc, char **argv) {
    // Parse args
    args::ArgumentParser parser("Inverse convex hull filling");
    // make parser for input name and other params
    args::HelpFlag help(parser,    "help", "Display this help menu", {'h', "help"});
    // args::Flag do_just_ours(parser, "do_just_ours", "do just ours", {"just_ours"});
    // args::ValueFlag<int> total_samples(parser, "ICOS_samples", "Total number of samples", {"samples"});
    args::ValueFlag<std::string> convex_input_arg(  parser, "convex_input_str", "path to input convex shape to be filled", {'c', "convex_dir"});
    args::ValueFlag<std::string> concave_input_arg( parser, "concave_input_arg", "path to concave shape to fill with", {'m', "concave_dir"});
    args::ValueFlag<std::string> remesh_scale_arg( parser, "remesh_scale_arg", "Remesh scale arg", {"remesh"});
    args::ValueFlag<std::string> deformed_input_arg(parser, "deformed_input_arg", "path to deformed shape that fills hull", {"deformed_dir"});
    args::ValueFlag<std::string> goal_G_text_arg(   parser, "goal_G_text_arg", "path to text containing G", {'g', "com"});

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

    if (remesh_scale_arg){
        concave_input_remesh_scale = std::stod(args::get(remesh_scale_arg));
    }

    std::cout << ANSI_FG_YELLOW << "input convex shape: " << hull_input_name << ANSI_RESET << "\n";
    std::cout << ANSI_FG_YELLOW << "input concave shape: " << concave_input_name << ANSI_RESET << "\n";
    std::cout << ANSI_FG_YELLOW << "goal G: " << goal_G << ANSI_RESET << "\n";

    // build mesh
    vis_utils = VisualUtils();
    polyscope::init();
    if (!deformed_input_arg){
        initialize_state(hull_input_name, concave_input_name, true);
    }
    else{
        initialize_state(hull_input_name, concave_input_name, false);
        // polyscope::show();
        std::string deformed_input_path = args::get(deformed_input_arg);
        std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
        std::unique_ptr<VertexPositionGeometry> geometry_ptr;
        std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(deformed_input_path);
        deformed_mesh = mesh_ptr.release();
        deformed_geometry = geometry_ptr.release();
        // TODO make sure input concave geo comes from the same mesh
        if (concave_mesh->nVertices() != deformed_mesh->nVertices()){
            std::cerr << "INPUT: deformed mesh and concave mesh have different number of vertices\n";
            return 1;
        }

        // visuals
        polyscope::registerSurfaceMesh("deformed input",
                                    deformed_geometry->inputVertexPositions, deformed_mesh->getFaceVertexList())->setTransparency(0.7);
        // check probabilities for the deformed shape
        Forward3DSolver* deformed_solver = new Forward3DSolver(deformed_mesh, deformed_geometry, goal_G, true); //
        deformed_solver->set_uniform_G();
        deformed_solver->initialize_pre_computes();
        BoundaryBuilder* deformed_bnd_builder = new BoundaryBuilder(deformed_solver);
        deformed_bnd_builder->build_boundary_normals();
        std::cout << "deformed G   :" << deformed_solver->get_G() << "\n";
        std::cout << "goal G was   : " << goal_G << "\n";

        // visuals
        // visualize_gauss_map(deformed_solver);
        // vis_utils.update_visuals(deformed_solver, deformed_bnd_builder, sphere_mesh, sphere_geometry);
        // vis_utils.draw_stable_patches_on_gauss_map(false, deformed_bnd_builder, false);
        // polyscope::show();
    }

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
            DeformationSolver* deformationSolver = 
                        new DeformationSolver(concave_mesh, concave_geometry, hull_mesh, hull_geometry, goal_G);   
            initialize_deformation_params(deformationSolver);
            size_t current_fill_iter = 0;
            Eigen::MatrixXd new_points = deformationSolver->solve_for_bending(1);

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
            // None yet

            // ideal stats
            final_solver->set_G(goal_G);
            final_solver->initialize_pre_computes();
            BoundaryBuilder* ideal_bnd_builder = new BoundaryBuilder(final_solver);
            ideal_bnd_builder->build_boundary_normals();
            std::cout << "post-deform with ideal G probabilities:\n";
            ideal_bnd_builder->print_area_of_boundary_loops();
            // visuals for post deform
            // visuals
            visualize_gauss_map(final_solver);
            vis_utils.update_visuals(final_solver, ideal_bnd_builder, sphere_mesh, sphere_geometry);
            vis_utils.draw_stable_patches_on_gauss_map(false, ideal_bnd_builder, false);
        }
        else if (animate_G_deform){
            polyscope::getPointCloud("goal G")->setEnabled(false);
            polyscope::getPointCloud("current G")->setEnabled(false);
            animate_G_deform = false;
            DeformationSolver* deformationSolver = 
                        new DeformationSolver(concave_mesh, concave_geometry, deformed_geometry, hull_mesh, hull_geometry, goal_G);   
            initialize_deformation_params(deformationSolver);
            Eigen::MatrixXd new_points = deformationSolver->solve_for_G(1);

            optimized_concave_mesh = deformationSolver->mesh;
            optimized_concave_geometry = deformationSolver->deformed_geometry;
            
            // checking probabilities for the deformed shape
            Forward3DSolver* final_solver = new Forward3DSolver(deformationSolver->mesh, deformationSolver->deformed_geometry, goal_G, true); // 
            
            polyscope::registerSurfaceMesh("deformed final hull", final_solver->hullGeometry->inputVertexPositions, 
                                           final_solver->hullMesh->getFaceVertexList())->setTransparency(0.4);

            
            // post process
            // ideal G probs
            final_solver->set_G(goal_G);
            final_solver->initialize_pre_computes();
            BoundaryBuilder* tmp_bnd_builder = new BoundaryBuilder(final_solver);
            tmp_bnd_builder->build_boundary_normals();
            std::cout << "  **  assuming perfect G probs:\n";
            tmp_bnd_builder->print_area_of_boundary_loops(); 

            // probs of new deformed shape
            final_solver->set_uniform_G();
            Vector3 post_deform_G = final_solver->get_G();
            std::cout << "goal G was   : " << goal_G << "\n";
            std::cout << "post-deform G: " << post_deform_G << "\n";
            
            final_solver->initialize_pre_computes();
            tmp_bnd_builder = new BoundaryBuilder(final_solver);
            tmp_bnd_builder->build_boundary_normals();

            
            printf(" post deform probabilities:\n");
            tmp_bnd_builder->print_area_of_boundary_loops();

            // visuals
            visualize_gauss_map(final_solver);
            vis_utils.update_visuals(final_solver, tmp_bnd_builder, sphere_mesh, sphere_geometry);
            vis_utils.draw_stable_patches_on_gauss_map(false, tmp_bnd_builder, false);
        }
        polyscope::frameTick();
    }
    polyscope::show();

    return EXIT_SUCCESS;
}


//