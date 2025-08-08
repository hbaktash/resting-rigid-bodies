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
// PhysicsEnv* my_env;
double bullet_step_size = 0.01; // 0.016;

int step_count = 1;
Vector3 G;
Vector3 pre_shift_G;

float orientation_gui[3] = {0., -1., 0.};

// stuff for Gauss map
float face_normal_vertex_gm_radi = 0.03,
      gm_distance = 2.,
      gm_radi = 1.;
int arcs_seg_count = 13;
Vector3 shift = {0., gm_distance , 0.},
        colored_shift = {gm_distance, gm_distance , 0.};
float arc_curve_radi = 0.01;
double friction_coeff = 0.1;

double ground_box_y = -2.1;
Vector3 ground_box_shape({10,1,10});

// Vector3 default_face_color({0.99,0.99,0.99});
Vector3 default_face_color({240./256.,178/256.,44./256.});

//
float scale_for_save = 1.;

// example choice
std::vector<std::string> all_polyhedra_items = {std::string("tet"), std::string("tet2"), std::string("tet0"),std::string("cube"), std::string("tilted cube"), std::string("sliced tet"), std::string("worst_case"), std::string("fox"), std::string("small_bunny"), std::string("bunnylp"), std::string("kitten"), std::string("double-torus"), std::string("knuckle_bone_real"),std::string("soccerball"), std::string("bunny"), std::string("gomboc"), std::string("dragon1"), std::string("dragon3"), std::string("mark_gomboc"), std::string("KnuckleboneDice"), std::string("Duende"), std::string("papa_noel"), std::string("reno"), std::string("baby_car"), std::string("rubberDuckie")};
std::string all_polygons_current_item = "bunny";

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


//GM stuff
ManifoldSurfaceMesh* sphere_mesh;
VertexPositionGeometry* sphere_geometry;
bool gm_is_drawn = false;
bool draw_snail_trail = true,
     animate = true;
// bool save_pos_to_file = false;
// Vector3 old_g_vec, new_g_vec;
int snail_trail_dummy_counter = 0;

int solid_angle_face_index = 0;


// for file saving
std::string mesh_title = "mesh";
bool save_original = false,
     save_oriented = false,
     save_orientation_trail_to_file = false,
     save_orientation_trail_scrs = false,
     snapshot_trail = false,
     save_quasi_trail_to_file = false,
     snapshot_quasi_trail = false;

int stable_face_index = 0; // index of the stable face to visualize
bool save_stable_orientation_pose = false;

std::string multi_orientations_path;
// Functions 


void initialize_vis(bool with_plane = true){
    polyscope::registerSurfaceMesh("my polyhedra", geometry->inputVertexPositions, mesh->getFaceVertexList());
    // ground plane on Polyscope has a weird height setting (scaled factor..)
    if (with_plane){
      auto psPlane = polyscope::addSceneSlicePlane("ground plane");
      psPlane->setDrawPlane(true);  // render the semi-transparent gridded plane
      psPlane->setDrawWidget(false);
      psPlane->setPose(glm::vec3{0., ground_box_y + 1, 0.}, glm::vec3{0., 1., 0.});
    }
}


Vector3 orientation_from_string(std::string filename){
    // the last part of the path is the orientation
        std::string orientation_str = filename.substr(filename.find_last_of("/\\") + 1);
        // there is no extension
        double x, y, z;
        // now split by '_'; just do substring 3 times
        size_t pos1 = orientation_str.find('_');
        size_t pos2 = orientation_str.find('_', pos1 + 1);
        x = std::stod(orientation_str.substr(0, pos1));
        y = std::stod(orientation_str.substr(pos1 + 1, pos2 - pos1 - 1));
        z = std::stod(orientation_str.substr(pos2 + 1));
        std::cout << "orientation from file: " << x << " " << y << " " << z << "\n";
        return Vector3({x, y, z});
  }


void generate_polyhedron_example(std::string mesh_full_path, bool preprocess = true){
  
  std::unique_ptr<SurfaceMesh> nm_mesh_ptr;
  std::unique_ptr<VertexPositionGeometry> nm_geometry_ptr;
  std::tie(nm_mesh_ptr, nm_geometry_ptr) = readSurfaceMesh(mesh_full_path);
  SurfaceMesh *nm_mesh = nm_mesh_ptr.release();
  VertexPositionGeometry *nm_geometry = nm_geometry_ptr.release();
  nm_mesh->greedilyOrientFaces();
  // for (Vertex v: nm_mesh->vertices()){
  //   if (!v.isManifold()){
  //     std::cout << "non-manifold vertex "<< v.getIndex() << "\n";
  //   }
  // }
  // for (Edge e: nm_mesh->edges()){
  //   if (!e.isManifold()){
  //     std::cout << "non-manifold edge "<< e.getIndex()<< ": " << e.firstVertex().getIndex() << " " << e.secondVertex().getIndex() << "\n";
  //     std::cout << "sibling: " << e.halfedge().sibling().getIndex() << ": "
  //               << e.halfedge().sibling().tailVertex().getIndex() << " " 
  //               << e.halfedge().sibling().tipVertex().getIndex() << "\n";
  //     std::cout << "sibling.sibling: " << e.halfedge().sibling().sibling().getIndex() << ": "
  //               << e.halfedge().sibling().sibling().tailVertex().getIndex() << " " 
  //               << e.halfedge().sibling().sibling().tipVertex().getIndex() << "\n";
  //   }
  // }
  nm_mesh->compress();

  mesh_ptr = nm_mesh->toManifoldMesh();
  mesh = mesh_ptr.release();
  geometry = new VertexPositionGeometry(*mesh);
  // transfer from nm geometry
  for (Vertex v : mesh->vertices()) {
    geometry->inputVertexPositions[v.getIndex()] = nm_geometry->inputVertexPositions[v.getIndex()];
  }

  // initial rotation for piggy
  for (Vertex v: mesh->vertices()){
    Vector3 p = geometry->inputVertexPositions[v];
    geometry->inputVertexPositions[v] = Vector3({-p.x, p.z, p.y});
  }

  // preproccess and shift for external use
  if (preprocess){
    preprocess_mesh(mesh, geometry, true, false);
    G = find_center_of_mass(*mesh, *geometry).first;
    pre_shift_G = G;
    // std::cout << "center of mass before shift: " << G << "\n";
    for (Vertex v: mesh->vertices()){
      geometry->inputVertexPositions[v] -= G;
    }
    G = find_center_of_mass(*mesh, *geometry).first;
  }
  G = find_center_of_mass(*mesh, *geometry).first;
  // std::cout << "center of mass after shift: " << G << "\n";
  // double max_dist = 0;
  // for (Vertex v: mesh->vertices()){
  //   max_dist = std::max(max_dist, geometry->inputVertexPositions[v].norm());
  // }
  // std::cout << "max dist from center: " << max_dist << "\n";
}


void update_solver(){
  forwardSolver = new Forward3DSolver(mesh, geometry, G, true);
  forwardSolver->initialize_pre_computes();
  boundary_builder = new BoundaryBuilder(forwardSolver);
  boundary_builder->build_boundary_normals();
}

void color_faces(Forward3DSolver *fwd_solver){
  // printf("hull faces: %d\n", forwardSolver->hullMesh->nFaces());
  face_colors = FaceData<Vector3>(*fwd_solver->hullMesh, default_face_color);
  std::vector<Face> stable_faces;
  for (Face f: fwd_solver->hullMesh->faces()){
    if (fwd_solver->face_is_stable(f))
      stable_faces.push_back(f);
  }
  face_colors = generate_random_colors(fwd_solver->hullMesh, stable_faces);
  
  for (Face f: fwd_solver->hullMesh->faces()){
    if (!fwd_solver->face_is_stable(f))
      face_colors[f] = default_face_color;
  }
}


void visualize_gauss_map(Forward3DSolver* forwardSolver){
    std::vector<Vector3> normals_icosahedral = generate_normals_icosahedral(40);
    // get the convex hull of the normals
    std::tie(sphere_mesh, sphere_geometry) = get_convex_hull_mesh(normals_icosahedral);
    //update vis utils
    vis_utils.draw_gauss_map(forwardSolver, sphere_mesh, sphere_geometry);
}


void init_visuals(){
  auto psInputMesh = polyscope::registerSurfaceMesh(
    "init input mesh",
    geometry->inputVertexPositions, mesh->getFaceVertexList());
  psInputMesh->setTransparency(1.0);
  psInputMesh->setEnabled(true);
  auto psHullMesh = polyscope::registerSurfaceMesh(
    "init hull mesh",
    forwardSolver->hullGeometry->inputVertexPositions, forwardSolver->hullMesh->getFaceVertexList(),
    polyscopePermutations(*forwardSolver->hullMesh));
  psHullMesh->setEnabled(true);
  psHullMesh->setEdgeWidth(0.7);
  psHullMesh->setTransparency(0.9);
  vis_utils.draw_G(forwardSolver->get_G());
  visualize_gauss_map(forwardSolver);
}

void update_solver_and_boundaries(){
  auto t1 = clock();
  forwardSolver->set_G(G);
  // printf("forward precomputes \n");
  forwardSolver->initialize_pre_computes();
  printf("building boundary normals \n");
  boundary_builder->build_boundary_normals();
}


void draw_stable_patches_on_gauss_map(bool on_height_surface = false){
  forwardSolver->initialize_pre_computes();
  // std::vector<Vector3> boundary_normals;
  auto net_pair = build_and_draw_stable_patches_on_gauss_map(boundary_builder,
                                                              vis_utils.center, vis_utils.gm_radi, vis_utils.arcs_seg_count, 
                                                              on_height_surface);
}

void draw_trail_on_gm(std::vector<Vector3> trail, glm::vec3 color, std::string name, double radi, bool color_gradient = false){
  std::vector<std::pair<size_t, size_t>> edge_inds;
  if (!color_gradient){
    for (size_t i = 0; i < trail.size()-1; i++)
      edge_inds.push_back({i, i+1});
    draw_arc_network_on_sphere(edge_inds, trail, vis_utils.center, vis_utils.gm_radi, vis_utils.arcs_seg_count, name, radi, color);//{0.7,0.1,0.8}
  }
  else{
    std::vector<glm::vec3> colors;
    for (size_t i = 0; i < trail.size(); i++){
      glm::vec3 tmp_color = color;
      tmp_color.x = tmp_color.x * ((i/(float)trail.size()));
      draw_arc_on_sphere(trail[i], trail[i+1], vis_utils.center, vis_utils.gm_radi, vis_utils.arcs_seg_count, i, radi, tmp_color);
    }
  }
}


void initialize_env(Vector3 orientation, bool visuals = true){
    // physics env
    my_env = new PhysicsEnv();
    my_env->init_physics();
    my_env->init_geometry(forwardSolver->hullMesh, forwardSolver->hullGeometry);
    my_env->add_ground(ground_box_y, ground_box_shape);
    my_env->add_object(G, orientation);
    my_env->default_step_size = bullet_step_size;
    // // polyscope
    // if (visuals)
    //     initialize_vis(true);
}

void build_quasi_static_snail_trail(Vector3 initial_orientation, glm::vec3 color, std::string name, double radi = 0.01){
    std::vector<Vector3> snail_trail = forwardSolver->snail_trail_log(initial_orientation);
    // S^2 shift
    draw_trail_on_gm(snail_trail, color, name, radi);
}

std::pair<std::vector<Eigen::Matrix4d>, std::vector<Vector3>> 
generate_transformations_for_orientation_sequence(
    Vector3 initial_orientation,
    Forward3DSolver* forwardSolver,
    Vector3 floor_vec,
    double goal_angle_step = 0.005
){
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
quasi_static_snail_trail_to_global_transformations(Vector3 initial_orientation, bool visualize = false){
    // initialize the trail in forwardSolver
    std::vector<Vector3> snail_trail = forwardSolver->snail_trail_log(initial_orientation);
    // split the snail trail
    auto [local_transformation_matrices, saved_snail_trail_refined] = generate_transformations_for_orientation_sequence(
        initial_orientation, forwardSolver, Vector3({0,-1,0}), 0.005);
    // local to global
    std::vector<Eigen::Matrix4d> transformation_matrices;
    Eigen::Matrix4d global_transformation = Eigen::Matrix4d::Identity();
    for (int i = 0; i < local_transformation_matrices.size(); i++){
        Eigen::Matrix4d local_transformation = local_transformation_matrices[i];
        global_transformation = local_transformation * global_transformation;
        transformation_matrices.push_back(global_transformation);
    }
    // apply to the mesh and visualize each step
    if (visualize){
        std::cout << "visualizing snail trail steps \n";
        VertexData<Vector3> tmp_positions(*forwardSolver->inputMesh);
        // tmp_positions = forwardSolver->inputGeometry->inputVertexPositions;
        Vector3 tmp_orientation = initial_orientation.normalize();
        for (int i = 0; i < transformation_matrices.size(); i++){
            Eigen::Matrix4d transformation_matrix = transformation_matrices[i];
            // apply to the input mesh
            for (Vertex v: forwardSolver->inputMesh->vertices()){
                Vector3 p = forwardSolver->inputGeometry->inputVertexPositions[v];
                // apply the transformation matrix
                Eigen::Vector4d tmp_v({p.x,
                                       p.y,
                                       p.z, 
                                       1.0});
                tmp_v = transformation_matrix * tmp_v;
                tmp_positions[v] = Vector3({tmp_v[0], tmp_v[1], tmp_v[2]});
            }
            // visualize
            polyscope::registerSurfaceMesh("quasi static snail step ", 
                tmp_positions, 
                forwardSolver->inputMesh->getFaceVertexList())->setEnabled(true);
            // only the snail trail
            tmp_orientation = saved_snail_trail_refined[i];
            std::vector<Vector3> tmp_orientation_vec = {tmp_orientation * vis_utils.gm_radi + vis_utils.center};
            polyscope::registerPointCloud("quasi orientation on gm", 
                                            tmp_orientation_vec)->setPointRadius(0.03, false)->setPointColor({0,0,0});
            polyscope::show();
        }
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
    if (ImGui::SliderFloat3("orientation", orientation_gui, -10.0f, 10.0f)){
        Vector3 orientation({orientation_gui[0], orientation_gui[1], orientation_gui[2]});
        std::vector<Vector3> orientation_vec = {orientation.normalize() * vis_utils.gm_radi + vis_utils.center};
        polyscope::registerPointCloud("orientation on gm", orientation_vec)->setPointRadius(0.03, false);
    }
    ImGui::Checkbox("save original", &save_original);
    ImGui::Checkbox("save oriented", &save_oriented);
    if (ImGui::Button("visualize orientation ")) {
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
        
        
        // save to file
        // save the rotated mesh, orientation and center of mass after rotation to files
        std::string PARENT_SAVE_DIR = "/Users/hbakt/Library/CloudStorage/Box-Box/Rolling-Dragon presentations/Talk/Media/orientation_on_S2_slide";
        if (save_original){
            // save the original mesh first
            writeSurfaceMesh(*mesh, *geometry, PARENT_SAVE_DIR + "/" + mesh_title + "/" + mesh_title + "_original.obj");
        }
        if (save_oriented){
            // save oriented
            VertexPositionGeometry geo_for_save(*mesh);
            for (Vertex v: mesh->vertices()){
                geo_for_save.inputVertexPositions[v] = new_positions[v];
            }
            // convert orientation to string
            std::string orientation_str = std::to_string(orientation.x) + "_" + std::to_string(orientation.y) + "_" + std::to_string(orientation.z);
            std::string path_to_dir = PARENT_SAVE_DIR + "/" + mesh_title + "/orientations/" + orientation_str;
            // create directory if it does not exist
            if (!fs::exists(path_to_dir)){
                fs::create_directories(path_to_dir);
            }
            writeSurfaceMesh(*mesh, geo_for_save, path_to_dir + "/" + "oriented_shape.obj");
        }
    }
    ImGui::Checkbox("save orientation trail", &save_orientation_trail_to_file);
    ImGui::Checkbox("snapshot trail", &snapshot_trail);
    if (ImGui::Button("bullet snail trail")) {
        // initialize bullet env
        Vector3 orientation({orientation_gui[0], orientation_gui[1], orientation_gui[2]});
        orientation = orientation.normalize();
        initialize_env(orientation, false);
        // visualize snail trail 
        Face touching_face = my_env->final_stable_face(draw_snail_trail);
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
        
        if (save_orientation_trail_to_file){
            // save orientation trajectory to file
            std::string PARENT_SAVE_DIR = "/Users/hbakt/Library/CloudStorage/Box-Box/Rolling-Dragon presentations/Talk/Media/orientation_on_S2_slide";
            std::string orientation_str = std::to_string(orientation.x) + "_" + std::to_string(orientation.y) + "_" + std::to_string(orientation.z);
            std::string path_to_dir = PARENT_SAVE_DIR + "/" + mesh_title + "/orientations/" + orientation_str;
            // create directory if it does not exist
            if (!fs::exists(path_to_dir)){
                fs::create_directories(path_to_dir);
            }
            // save
            // save orientation trail
            std::ofstream ofs(path_to_dir + "/" + "orientation_trail.txt");
            for (const auto& pos : my_env->orientation_trail) {
                ofs << pos.x << " " << pos.y << " " << pos.z << "\n";
            }
            ofs.close();
            // save transformation matrices
            std::ofstream ofs_mat(path_to_dir + "/" + "transformation_matrices.txt");
            for (const auto& mat : my_env->trans_mat_trail) {
                ofs_mat << mat << "\n";
            }
            ofs_mat.close();
            // save height log
            std::ofstream ofs_height(path_to_dir + "/" + "height_log.txt");
            for (const auto& height : height_log) {
                ofs_height << height << "\n";
            }
            ofs_height.close();
        }

    }
    ImGui::Checkbox("save quasi static trail to file", &save_quasi_trail_to_file);
    if (ImGui::Button("Build quasistatic snail trail")){
        // get rotations corresponding to the snail trail
        Vector3 orientation({orientation_gui[0], orientation_gui[1], orientation_gui[2]});
        orientation = orientation.normalize();
        
        build_quasi_static_snail_trail(orientation, {39./255., 189./255., 0}, "quasi-static trail", 3.);
        auto [transformation_matrices, saved_snail_trail_refined] = quasi_static_snail_trail_to_global_transformations(orientation, false);
        
        // save to file like bullet snail trail matrices
        std::string PARENT_SAVE_DIR = "/Users/hbakt/Library/CloudStorage/Box-Box/Rolling-Dragon presentations/Talk/Media/orientation_on_S2_slide";
        std::string orientation_str = std::to_string(orientation.x) + "_" + std::to_string(orientation.y) + "_" + std::to_string(orientation.z);
        std::string path_to_dir = PARENT_SAVE_DIR + "/" + mesh_title + "/orientations/" + orientation_str + "/quasi_static_trail";
        if (save_quasi_trail_to_file){
            // create directory if it does not exist
            if (!fs::exists(path_to_dir)){
                fs::create_directories(path_to_dir);
            }
            // save transformation matrices
            std::ofstream ofs_mat(path_to_dir + "/" + "transformation_matrices.txt");
            for (const auto& mat : transformation_matrices) {
                ofs_mat << mat << "\n";
            }
            ofs_mat.close();
        }
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
    args::ValueFlag<std::string> orientation_path_arg(parser, "orientation", "path to orientation files", {"orientation"});


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
    std::string mesh_path, com_path, orientation_path;
    if (mesh_path_arg)
        mesh_path = args::get(mesh_path_arg);
    if (center_of_mass_path_arg)
        com_path = args::get(center_of_mass_path_arg);
    if (orientation_path_arg){
        orientation_path = args::get(orientation_path_arg);
        Vector3 orientation_inp = orientation_from_string(orientation_path);
        orientation_gui[0] = orientation_inp.x;
        orientation_gui[1] = orientation_inp.y;
        orientation_gui[2] = orientation_inp.z;
    }

    // extract mesh title
    mesh_title = mesh_path.substr(mesh_path.find_last_of("/\\") + 1);
    mesh_title = mesh_title.substr(0, mesh_title.find_last_of("."));
    // --- Debugging --- //
    polyscope::init();
    vis_utils = VisualUtils();
    generate_polyhedron_example(mesh_path);
    // initial 
    update_solver();

    init_visuals();

    draw_stable_patches_on_gauss_map();

        

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
    polyscope::state::userCallback = myCallback;
    polyscope::show();
    return EXIT_SUCCESS;  
  }