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
#include "unsupported/Eigen/EulerAngles"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

// #include "LinearMath/btVector3.h"
// #include "btBulletDynamicsCommon.h"

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

#include "chrono"
// #include "ipc/ipc.hpp"

namespace chrono = std::chrono;
using clock_type = chrono::high_resolution_clock;
using seconds_fp = chrono::duration<double, chrono::seconds::period>;

using namespace geometrycentral;
using namespace geometrycentral::surface;

// simulation stuff
// PhysicsEnv* my_env;
float step_size = 0.001, // 0.016
      refresh_x, refresh_y, refresh_z,
      G_x, G_y, G_z;
int step_count = 1;
Vector3 G;
Vector3 refresh_orientation({0,-1,0});

// stuff for Gauss map
float face_normal_vertex_gm_radi = 0.03,
      gm_distance = 2.,
      gm_radi = 1.;
int arcs_seg_count = 13;
Vector3 shift = {0., gm_distance , 0.},
        colored_shift = {gm_distance, gm_distance , 0.};
float arc_curve_radi = 0.01;


double ground_box_y = -2.1;
Vector3 ground_box_shape({10,1,10});

// Vector3 default_face_color({0.99,0.99,0.99});
Vector3 default_face_color({240./256.,178/256.,44./256.});

// example choice
std::vector<std::string> all_polyhedra_items = {std::string("tet"), std::string("tet2"), std::string("tet0"),std::string("cube"), std::string("tilted cube"), std::string("sliced tet"), std::string("worst_case"), std::string("fox"), std::string("small_bunny"), std::string("bunnylp"), std::string("kitten"), std::string("double-torus"), std::string("knuckle_bone_real"),std::string("soccerball"), std::string("bunny"), std::string("gomboc"), std::string("dragon1"), std::string("dragon3"), std::string("mark_gomboc"), std::string("KnuckleboneDice"), std::string("Duende"), std::string("papa_noel"), std::string("reno"), std::string("baby_car"), std::string("rubberDuckie")};
std::string all_polygons_current_item = "bunny";

// GC stuff
std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr;
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;


// raster image stuff
FaceData<Vector3> face_colors;
int total_samples = 1e4, 
    max_steps_IPC = 1000;
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

void generate_polyhedron_example(std::string poly_str, bool triangulate = false){
    // readManifoldSurfaceMesh()
  
  // if(poly_str == "baby_car"){
  //   printf("reading baby car\n");
  if (poly_str == "tet" || poly_str == "tet2" || poly_str == "sliced tet"){
    std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(poly_str);
    mesh = mesh_ptr.release();
    geometry = geometry_ptr.release();
    preprocess_mesh(mesh, geometry, triangulate || std::strcmp(poly_str.c_str(), "gomboc") == 0, false, 1.);
  }
  else{
    std::unique_ptr<SurfaceMesh> nm_mesh_ptr;
    std::unique_ptr<VertexPositionGeometry> nm_geometry_ptr;
    std::string filename = "../meshes/"+poly_str+".obj";
    std::tie(nm_mesh_ptr, geometry_ptr) = readSurfaceMesh(filename);
    SurfaceMesh *nm_mesh = nm_mesh_ptr.release();
    mesh_ptr = nm_mesh->toManifoldMesh();
    mesh = mesh_ptr.release();
    geometry = geometry_ptr.release();
    preprocess_mesh(mesh, geometry, true);
  }
  G = find_center_of_mass(*mesh, *geometry).first;
  for (Vertex v: mesh->vertices()){
    geometry->inputVertexPositions[v] -= G;
  }
  G = find_center_of_mass(*mesh, *geometry).first;
  double max_dist = 0;
  for (Vertex v: mesh->vertices()){
    max_dist = std::max(max_dist, geometry->inputVertexPositions[v].norm());
  }
  std::cout << "center of mass after shift: " << G << "\n";
  std::cout << "max dist from center: " << max_dist << "\n";
  
  // write for IPC use
  writeSurfaceMesh(*mesh, *geometry, "../../rigid-ipc/meshes/centered_COMs/" + all_polygons_current_item + ".obj");
}



void initialize_boundary_builder(){
  boundary_builder = new BoundaryBuilder(forwardSolver);
  boundary_builder->build_boundary_normals();
  boundary_builder->print_area_of_boundary_loops();
}


void update_solver(){
  //assuming convex input here
  auto t1 = clock();
  forwardSolver = new Forward3DSolver(mesh, geometry, G, true);
  forwardSolver->initialize_pre_computes();
  printf("solver precompute time: %f\n", (clock() - t1)/(double)CLOCKS_PER_SEC);
  initialize_boundary_builder();
  vis_utils.forwardSolver = forwardSolver;
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

// another polyhedra for the sake of a good colored raster image
void visualize_colored_polyhedra(){
  VertexData<Vector3> shifted_positions(*forwardSolver->hullMesh);
  for (Vertex v: forwardSolver->hullMesh->vertices()){
    shifted_positions[v] = forwardSolver->hullGeometry->inputVertexPositions[v] + colored_shift;
  }
  auto coloredPsMesh = polyscope::registerSurfaceMesh("colored polyhedra", shifted_positions, forwardSolver->hullMesh->getFaceVertexList());
  // generate random colors and color the faces
  color_faces(forwardSolver);
  polyscope::SurfaceFaceColorQuantity *faceQnty = coloredPsMesh->addFaceColorQuantity("random face colors", face_colors);
  faceQnty->setEnabled(true);
  // add colors to the original polyhedra as well
  polyscope::SurfaceFaceColorQuantity *faceQnty2 = polyscope::getSurfaceMesh("init hull mesh")->addFaceColorQuantity("random face colors2", face_colors);
  faceQnty2->setEnabled(true);
}

void visualize_gauss_map(){
  // just draw the sphere next to the main surface
  std::unique_ptr<ManifoldSurfaceMesh> sphere_mesh_ptr;
  std::unique_ptr<VertexPositionGeometry> sphere_geometry_ptr;
  std::tie(sphere_mesh_ptr, sphere_geometry_ptr) = generate_polyhedra("sphere");
  sphere_mesh = sphere_mesh_ptr.release();
  sphere_geometry = sphere_geometry_ptr.release();
  
  //update vis utils
  vis_utils.sphere_geometry = sphere_geometry;
  vis_utils.sphere_mesh = sphere_mesh;

  vis_utils.draw_gauss_map();
  polyscope::getSurfaceMesh("height surface_func")->setEnabled(false);
  gm_is_drawn = true;
}


void init_visuals(){
  auto psInputMesh = polyscope::registerSurfaceMesh(
    "init input mesh",
    geometry->inputVertexPositions, mesh->getFaceVertexList(),
    polyscopePermutations(*mesh));
  psInputMesh->setTransparency(0.75);
  psInputMesh->setEnabled(true);
  auto psHullMesh = polyscope::registerSurfaceMesh(
    "init hull mesh",
    forwardSolver->hullGeometry->inputVertexPositions, forwardSolver->hullMesh->getFaceVertexList(),
    polyscopePermutations(*forwardSolver->hullMesh));
  psHullMesh->setEnabled(false);
  vis_utils.draw_G();
  visualize_colored_polyhedra();
  visualize_gauss_map();
  // psInputMesh->addFaceVectorQuantity("normals", geometry->faceNormals);
}

void update_solver_and_boundaries(){
  auto t1 = clock();
  forwardSolver->set_G(G);
  // printf("forward precomputes \n");
  forwardSolver->initialize_pre_computes();
  printf("building boundary normals \n");
  // boundary_builder->build_boundary_normals();
  boundary_builder->build_boundary_normals();
}




void draw_stable_patches_on_gauss_map(bool on_height_surface = false){
  forwardSolver->initialize_pre_computes();
  // std::vector<Vector3> boundary_normals;
  auto net_pair = build_and_draw_stable_patches_on_gauss_map(boundary_builder,
                                                              vis_utils.center, vis_utils.gm_radi, vis_utils.arcs_seg_count, 
                                                              on_height_surface);
}


// sample and raster with quasi static solver
void build_raster_image(){
  FaceData<std::vector<Vector3>> face_samples(*forwardSolver->hullMesh);
  std::vector<Vector3> final_orientations;
  std::vector<Vector3> falsh;
  int total_invalids = 0, total_samples = 0;
  auto t1 = clock();
  double avg_time = 0.;
  for (int i = 0; i < total_samples; i++){
    Vector3 random_orientation = {randomReal(-1,1), randomReal(-1,1), randomReal(-1,1)};
    if (random_orientation.norm() <= 1){
      auto t2 = clock();
      total_samples++;
      if (i % 2000 == 0)
        printf("$$$ at sample %d\n", i);
      random_orientation /= norm(random_orientation);
    //   my_env->refresh(G, random_orientation);
    //   Face touching_face = my_env->final_stable_face(); // Invalid if not close to a face normal; wtf?
      avg_time += (double)((clock() - t2)/(double)CLOCKS_PER_SEC);
    //   final_orientations.push_back(my_env->get_current_orientation() + vis_utils.center);
    //   if (touching_face.getIndex() == INVALID_IND){
    //     total_invalids++;
    //     falsh.push_back(random_orientation + vis_utils.center); //  shift for visualization
    //     continue;
    //   }
    //   face_samples[touching_face].push_back(random_orientation);
    }
  }
  printf(" --- sampling time: %f\n", (double)(clock() - t1));
  printf(" ### total invalid faces: %d/%d\n", total_invalids, total_samples);
  std::vector<Vector3> raster_positions,
                       raster_colors;
                       
  FaceData<double> accum_face_areas(*forwardSolver->hullMesh, 0.);
  printf("empirical probs:\n");
  size_t unstable_faces = 0;
  for (Face f: forwardSolver->hullMesh->faces()){
    if (forwardSolver->face_last_face[f]!= f)
      unstable_faces++;
    std::vector<Vector3> tmp_points = face_samples[f];
    double prob = tmp_points.size()/(double)(total_samples - total_invalids);
    accum_face_areas[forwardSolver->face_last_face[f]] += prob;
    if(tmp_points.size() != 0) printf(" --- f %d: %f\n", f.getIndex(), prob);
    for (Vector3 tmp_p: tmp_points){
      raster_positions.push_back(tmp_p + shift);
      raster_colors.push_back(face_colors[forwardSolver->face_last_face[f]]); // 
      // std::cout<< "tmp color is: " << face_colors[f] << "\n";
    }
  }
  printf(" &&& unstable faces: %d\n", unstable_faces);
  printf("accumulated empirical -VS- MS complex:\n");
  for (Face f: forwardSolver->hullMesh->faces()){
    if (forwardSolver->face_is_stable(f)){
      printf("  f %d -> emp: %f    ----   MS: %f \n", f.getIndex(), accum_face_areas[f], boundary_builder->face_region_area[f]/(4.*PI));
    }
  }
  //
  printf(" $$$ KL divergence! $$$\n");
  double kl_divergence = 0.;
  for (Face f: forwardSolver->hullMesh->faces()){
    if (forwardSolver->face_is_stable(f)){
      double emp_prob = accum_face_areas[f];
      double ms_prob = boundary_builder->face_region_area[f]/(4.*PI);
      if (emp_prob != 0. && ms_prob != 0.)
        kl_divergence += ms_prob * std::log(ms_prob/emp_prob);
    }
  }
  printf("KL divergence: %f\n", kl_divergence);

  auto final_ori_pc = polyscope::registerPointCloud("final orientations", final_orientations);
  final_ori_pc->setPointColor({0.1,0.1,0.1});
  final_ori_pc->setPointRadius(0.04, false);
  final_ori_pc->setEnabled(true);

  auto raster_pc = polyscope::registerPointCloud("raster point cloud", raster_positions);
  polyscope::PointCloudColorQuantity* pc_col_quant = raster_pc->addColorQuantity("random color", raster_colors);
  pc_col_quant->setEnabled(true);
  auto falsh_pc = polyscope::registerPointCloud("falsh point cloud", falsh);
  falsh_pc->setPointColor({0.,0.,0.});
  falsh_pc->setEnabled(true);
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


Eigen::AngleAxisd aa_from_init_ori(Vector3 init_ori){
    Vector3 rotation_axis = cross(init_ori, Vector3({0,-1,0})).normalize();
    double rotation_angle = angle(init_ori, Vector3({0,-1,0})); // WARNING: uses acos

    Eigen::AngleAxisd Rinput_aa(rotation_angle, Eigen::Vector3d(rotation_axis.x, rotation_axis.y, rotation_axis.z));    
    return Rinput_aa;
}


Eigen::EulerAnglesZYXd euler_from_init_ori(Vector3 init_ori){
    Eigen::EulerAnglesZYXd temp_R_euler(aa_from_init_ori(init_ori).toRotationMatrix());
    return temp_R_euler;
}


Eigen::AngleAxisd run_sim_fetch_rot(Vector3 init_ori, nlohmann::json jf, int max_iters, 
                                    std::string exec_dir, std::string jsons_dir, std::string example_fname, std::string type){
  
  Eigen::AngleAxisd Rinput_aa = aa_from_init_ori(init_ori);
  Eigen::EulerAnglesZYXd temp_R_euler(Rinput_aa.toRotationMatrix());
  jf["max_iterations"] = max_iters;
  jf["rigid_body_problem"]["rigid_bodies"][0]["mesh"] = "centered_COMs/" + all_polygons_current_item + ".obj";
  jf["rigid_body_problem"]["rigid_bodies"][0]["rotation"] = {temp_R_euler.angles()[2] * 180. / M_PI, 
                                                              temp_R_euler.angles()[1] * 180. / M_PI, 
                                                              temp_R_euler.angles()[0] * 180. / M_PI};
  jf["rigid_body_problem"]["rigid_bodies"][0]["position"] = {0, 1, 0};

  std::string input_name = type + "_" + all_polygons_current_item + ".json";
  std::ofstream ofs(jsons_dir + "/" + input_name);
  ofs << jf;
  ofs.close();
  // auto t1 = clock_type::now();
  std::string cmd = exec_dir + 
                  " --chkpt 10001 --ngui --nthreads 1 " + // threads = 1 for the non-deteministic behavior
                  jsons_dir + "/" + input_name + " " + 
                  jsons_dir + "/temp_out";
  // std::cout << "\nrunning: ...\n" ;//<< cmd << "\n";
  int out = system((cmd + "> /dev/null").c_str()); //
  // std::cout << "result:" <<  out << "\n\n"; // suppress output
  // auto t2 = clock_type::now();
  // std::cout << "  &&& sample run time " << chrono::duration_cast<seconds_fp>(t2 - t1).count() << " seconds\n";
  
  // read output
  std::ifstream ifs_out(jsons_dir + "/temp_out/sim.json");
  nlohmann::json jf_out = nlohmann::json::parse(ifs_out);

  size_t num_states = jf_out["animation"]["state_sequence"].size();
  auto r0_aa = jf_out["animation"]["state_sequence"][0]["rigid_bodies"][0]["rotation"],
        rn_aa = jf_out["animation"]["state_sequence"][num_states - 1]["rigid_bodies"][0]["rotation"];
  Eigen::Vector3d r0_vec(r0_aa[0], r0_aa[1], r0_aa[2]),
                  rn_vec(rn_aa[0], rn_aa[1], rn_aa[2]);
  Eigen::AngleAxisd R0_aa(r0_vec.norm(), r0_vec.normalized()),
                    Rn_aa(rn_vec.norm(), rn_vec.normalized());
  Eigen::AngleAxisd R_inert_aa(Rinput_aa.inverse().toRotationMatrix() * R0_aa.toRotationMatrix());
  Eigen::AngleAxisd R_rest_aa(Rn_aa.toRotationMatrix() * R_inert_aa.inverse().toRotationMatrix()); 

  return R_rest_aa;
}


void run_IPC_samples_MCMC(std::string example_fname, int total_samples = 1, int max_iters = 200){
    std::string jsons_dir = "/Users/hbakt/Desktop/code/rigid-ipc/fixtures/3D/examples";
    std::string exec_dir = "/Users/hbakt/Desktop/code/rigid-ipc/build/rigid_ipc_sim";
    std::ifstream ifs(jsons_dir + "/" + example_fname);
    nlohmann::json jf = nlohmann::json::parse(ifs);
    ifs.close();

    Eigen::EulerAnglesXYZd R0_euler(0, 0, 0); // 0.5*M_PI
    FaceData<std::vector<Vector3>> face_to_ori(*forwardSolver->hullMesh);
    FaceData<size_t> face_counts(*forwardSolver->hullMesh);

    int samples = 0;
    size_t valid_count = 0;
    auto t1 = clock();
    while (samples < total_samples){
        Vector3 random_orientation = {randomReal(-1,1), randomReal(-1,1), randomReal(-1,1)};
        if (random_orientation.norm() <= 1){
            // sample random orientation -> euler rotation
            samples++;
            if (samples % 50 == 0){
                printf("$$$ at sample %d\n", samples);
                printf("avg time per sample: %f\n", (clock() - t1)/((double) samples * (double)CLOCKS_PER_SEC));          
                for (Face f: forwardSolver->hullMesh->faces()){
                    if (face_counts[f] != 0)
                        std::cout << "face " << f.getIndex() << " count: " << face_counts[f]/(double)valid_count << "\n";
                }
            }
            random_orientation = random_orientation.normalize();
            
            // run the sim
            Eigen::AngleAxisd R_rest_aa = run_sim_fetch_rot(random_orientation, jf, max_iters, exec_dir, jsons_dir, example_fname, "MC");
            VertexData<Vector3> rotated_poses(*forwardSolver->hullMesh);
            for (Vertex v: forwardSolver->hullMesh->vertices()){
                rotated_poses[v] = vec2vec3(R_rest_aa.toRotationMatrix() * vec32vec(forwardSolver->hullGeometry->inputVertexPositions[v]));
            }
            VertexPositionGeometry rotated_geo(*forwardSolver->hullMesh, rotated_poses);

            Face lowest_face;
            double lowest_height = 1e9;
            for (Face f: forwardSolver->hullMesh->faces()){
                Vector3 rotated_face_normal = rotated_geo.faceNormal(f);
                double height = cross(rotated_face_normal, Vector3({0,-1,0})).norm();
                if (height < lowest_height){
                    lowest_face   = f;
                    lowest_height = cross(rotated_face_normal, Vector3({0,-1,0})).norm();
                }
            }
            // for debugging
            Eigen::EulerAnglesZYXd temp_R_euler = euler_from_init_ori(random_orientation);
            if (lowest_height < 1e-3){
                valid_count++;
                face_counts[lowest_face]++;
                face_to_ori[lowest_face].push_back(random_orientation);
                
                if (forwardSolver->face_next_face[lowest_face] != lowest_face){
                    std::cout << "invalid orientation: " << temp_R_euler.angles().reverse().transpose() * 180. / M_PI << "\n";
                    std::cout << "with face: " << lowest_face.getIndex() << " ,  normal diff " << lowest_height << "\n";
                }
            }
            else{
                std::cout << "invalid orientation: " << temp_R_euler.angles().reverse().transpose() * 180. / M_PI << "\n";
            }
        }
    }
    printf("total time: %f\n", (double)(clock() - t1)/(double)CLOCKS_PER_SEC);
    printf("avg time per sample: %f\n", (clock() - t1)/((double) samples * (double)CLOCKS_PER_SEC));          
    for (Face f: forwardSolver->hullMesh->faces()){
        if (face_counts[f] != 0)
            std::cout << "face " << f.getIndex() << " count: " << face_counts[f]/(double)valid_count << "\n";
    }
    std::vector<Vector3> raster_positions,
                       raster_colors;
    for (Face f: forwardSolver->hullMesh->faces()){
        for (Vector3 tmp_p: face_to_ori[f]){
            raster_positions.push_back(tmp_p + shift);
            raster_colors.push_back(face_colors[f]);
        }
    }
    auto raster_pc = polyscope::registerPointCloud("raster pc IPC", raster_positions);
    polyscope::PointCloudColorQuantity* pc_col_quant = raster_pc->addColorQuantity("stable face colors", raster_colors);
    pc_col_quant->setEnabled(true);
    std::cout << "valid count: " << valid_count << "/" << total_samples << "\n";
    // json js = json::parse(str);
}


void run_IPC_samples_ICOS(std::string example_fname, int total_samples = 1, int max_iters = 200){
    std::string jsons_dir = "/Users/hbakt/Desktop/code/rigid-ipc/fixtures/3D/examples";
    std::string exec_dir = "/Users/hbakt/Desktop/code/rigid-ipc/build/rigid_ipc_sim";
    std::ifstream ifs(jsons_dir + "/" + example_fname);
    nlohmann::json jf = nlohmann::json::parse(ifs);
    ifs.close();
    Eigen::EulerAnglesXYZd R0_euler(0, 0, 0); // 0.5*M_PI
    FaceData<std::vector<Vector3>> face_to_ori(*forwardSolver->hullMesh);
    FaceData<size_t> face_counts(*forwardSolver->hullMesh);
    FaceData<double> face_dual_sum_areas(*forwardSolver->hullMesh);

    int samples = 0;
    int resolution = (int)sqrt(total_samples/10); // decent approximation
    size_t valid_count = 0;
    auto t1 = clock();
    ManifoldSurfaceMesh* icos_sphere_mesh;
    VertexPositionGeometry* icos_sphere_geometry;
    std::cout << "generating icosahedral sphere with resolution: " << resolution << "\n";
    std::tie(icos_sphere_mesh, icos_sphere_geometry) = get_convex_hull_mesh(generate_normals_icosahedral(resolution));

    // tilt the sphere randomly
    while (true){
        Vector3 random_orientation = {randomReal(-1,1), randomReal(-1,1), randomReal(-1,1)};
        if (random_orientation.norm() <= 1){
            random_orientation = random_orientation.normalize();
            Eigen::AngleAxisd R = aa_from_init_ori(random_orientation);
            for (Vertex v: icos_sphere_mesh->vertices())
                icos_sphere_geometry->inputVertexPositions[v] = vec2vec3(R.toRotationMatrix() * vec32vec(icos_sphere_geometry->inputVertexPositions[v]));
            break;
        }
    }
    
    // compute the total area of the ICOsphere
    double total_area = 0;
    for (Face f: icos_sphere_mesh->faces()){
        total_area += icos_sphere_geometry->faceArea(f);
    }
    // iterate over the sphere
    int verbose_period = 50;
    std::cout << "Icos subdivision vertex count N = " << icos_sphere_mesh->nVertices() << "\n";
    auto t0 = clock_type::now();
    auto first = clock_type::now();
    for (Vertex sample_v: icos_sphere_mesh->vertices()){
        Vector3 random_orientation = icos_sphere_geometry->inputVertexPositions[sample_v];
        double dual_area = icos_sphere_geometry->vertexDualArea(sample_v);
        random_orientation = random_orientation.normalize();
        
    
        // run the sim
        auto t1 = clock_type::now();
        Eigen::AngleAxisd R_rest_aa = run_sim_fetch_rot(random_orientation, jf, max_iters, exec_dir, jsons_dir, example_fname, "ICOS");
        // get the bottom face
        VertexData<Vector3> rotated_poses(*forwardSolver->hullMesh);
        for (Vertex v: forwardSolver->hullMesh->vertices()){
            rotated_poses[v] = vec2vec3(R_rest_aa.toRotationMatrix() * vec32vec(forwardSolver->hullGeometry->inputVertexPositions[v]));
        }
        VertexPositionGeometry rotated_geo(*forwardSolver->hullMesh, rotated_poses);
        std::cout << "      single instance time: " << chrono::duration_cast<seconds_fp>(clock_type::now() - t1).count() << " s\n";
          

        Face lowest_face;
        double lowest_height = 1e9;
        for (Face f: forwardSolver->hullMesh->faces()){
            Vector3 rotated_face_normal = rotated_geo.faceNormal(f);
            double height = cross(rotated_face_normal, Vector3({0,-1,0})).norm();
            if (height < lowest_height){
                lowest_face   = f;
                lowest_height = cross(rotated_face_normal, Vector3({0,-1,0})).norm();
            }
        }
        // for debugging
        Eigen::EulerAnglesZYXd temp_R_euler = euler_from_init_ori(random_orientation);    
        if (lowest_height < 1e-3){
            // std::cout << "stable face: " << f.getIndex() << "\n";
            valid_count++;
            face_counts[lowest_face]++;
            face_dual_sum_areas[lowest_face] += dual_area;
            face_to_ori[lowest_face].push_back(random_orientation);

            if (forwardSolver->face_next_face[lowest_face] != lowest_face){
                std::cout << "invalid orientation: " << temp_R_euler.angles().reverse().transpose() * 180. / M_PI << "\n";
                std::cout << "with face: " << lowest_face.getIndex() << " ,  normal diff " << lowest_height << "\n";
            }
        }
        else{
            std::cout << "invalid orientation: " << temp_R_euler.angles().reverse().transpose() * 180. / M_PI << "\n";
        }
        samples++;
        if (samples % verbose_period == 0){
          std::cout << "at sample: " << samples << "\n";
          auto last = clock_type::now();

          using seconds_fp = chrono::duration<double, chrono::seconds::period>;
          std::cout << "  -- average time: " << chrono::duration_cast<seconds_fp>(last - first).count()/(double)verbose_period << " seconds\n";
          first = last;
        }
    }
    std::cout << " ### TOTAL time: " << chrono::duration_cast<seconds_fp>(clock_type::now() - t0).count() << " seconds\n";
    std::cout << " ### average time per sample: " << chrono::duration_cast<seconds_fp>(clock_type::now() - t0).count()/(double)icos_sphere_mesh->nVertices() << " seconds\n";
    for (Face f: forwardSolver->hullMesh->faces()){
        if (face_counts[f] != 0)
            std::cout << "face " << f.getIndex() << " prob: " << face_dual_sum_areas[f]/total_area << "\n";
    }
    std::vector<Vector3> raster_positions,
                         raster_colors;
    for (Face f: forwardSolver->hullMesh->faces()){
        for (Vector3 tmp_p: face_to_ori[f]){
            raster_positions.push_back(tmp_p + shift);
            raster_colors.push_back(face_colors[f]);
        }
    }
    auto raster_pc = polyscope::registerPointCloud("raster pc IPC", raster_positions);
    polyscope::PointCloudColorQuantity* pc_col_quant = raster_pc->addColorQuantity("stable face colors", raster_colors);
    pc_col_quant->setEnabled(true);
    std::cout << "valid count: " << valid_count << "/" << icos_sphere_mesh->nVertices() << "\n";
    // json js = json::parse(str);
}

// to delete
void compare_quasi_sample_convergence(){
  
  // ICOS
  int resolution = (int)sqrt(total_samples/10);
  FaceData<double> ICOS_probs(*forwardSolver->hullMesh, 0.);
  int total_invalids = 0;
  size_t valid_count = 0;
  auto t1 = clock();
  ManifoldSurfaceMesh* icos_sphere_mesh;
  VertexPositionGeometry* icos_sphere_geometry;
  std::cout << "generating icosahedral sphere with resolution: " << resolution << "\n";
  std::tie(icos_sphere_mesh, icos_sphere_geometry) = get_convex_hull_mesh(generate_normals_icosahedral(resolution));
  // tilt the ICOS sphere randomly
  while (true){
      Vector3 random_orientation = {randomReal(-1,1), randomReal(-1,1), randomReal(-1,1)};
      if (random_orientation.norm() <= 1){
          random_orientation = random_orientation.normalize();
          Eigen::AngleAxisd R = aa_from_init_ori(random_orientation);
          for (Vertex v: icos_sphere_mesh->vertices())
              icos_sphere_geometry->inputVertexPositions[v] = vec2vec3(R.toRotationMatrix() * vec32vec(icos_sphere_geometry->inputVertexPositions[v]));
          break;
      }
  }
  // iterate over the ICOS sphere
  double total_area = 0;
  for (Face f: icos_sphere_mesh->faces()) total_area += icos_sphere_geometry->faceArea(f);
  std::cout << "start sampling... N = " << icos_sphere_mesh->nVertices() << "\n";
  total_samples = icos_sphere_mesh->nVertices();
  for (Vertex sample_v: icos_sphere_mesh->vertices()){
      Vector3 random_orientation = icos_sphere_geometry->inputVertexPositions[sample_v];
      double dual_area = icos_sphere_geometry->vertexDualArea(sample_v);
      random_orientation = random_orientation.normalize(); // redundant
      
      // run the quasi dyno
      Face touching_face = forwardSolver->final_touching_face(random_orientation);
      if (touching_face.getIndex() == INVALID_IND){
        total_invalids++;
        printf(" ----- WTFF invalid ICOS ----- !!\n");
        continue;
      }
      ICOS_probs[touching_face] += dual_area;
  }

  // ----------  Uniform rejection ----------  //

  FaceData<double> uniform_probs(*forwardSolver->hullMesh, 0.);
  int tmp_count = 0;
  total_invalids = 0;
  while (true){
    Vector3 random_g_vec = {randomReal(-1,1), randomReal(-1,1), randomReal(-1,1)};
    if (random_g_vec.norm() <= 1.){
      tmp_count++;
      if (tmp_count >= total_samples + 1)
        break;
      random_g_vec /= norm(random_g_vec);
      Face touching_face = forwardSolver->final_touching_face(random_g_vec);
      if (touching_face.getIndex() == INVALID_IND){
        total_invalids++;
        printf(" ----- WTFF invalid uniform ----- !!\n");
        continue;
      }
      uniform_probs[touching_face] += 1.;
    }
  }
  printf(" total samples uniform: %d\n", tmp_count);
  printf(" ICOS total samples uniform: %d\n", total_samples);
  double uni_diff_squared = 0., ICOS_diff_squared = 0.;
  for (Face f: forwardSolver->hullMesh->faces()){
    if (uniform_probs[f] != 0. || ICOS_probs[f] != 0.){
      double uniform_prob = uniform_probs[f] / (double)total_samples;
      double ICOS_prob = ICOS_probs[f]/total_area;
      double analytic_prob = boundary_builder->face_region_area[f]/(4.*PI);
      printf("f %d  anl: %f\n", f.getIndex(), analytic_prob);
      printf("      uniform: %f\n", uniform_prob);
      printf("          dif: %f\n", abs(uniform_prob - analytic_prob));
      printf("      ICOS: %f\n", ICOS_prob);
      printf("          dif: %f\n", abs(ICOS_prob - analytic_prob));
      uni_diff_squared += (uniform_prob - analytic_prob) * (uniform_prob - analytic_prob);
      ICOS_diff_squared += (ICOS_prob - analytic_prob) * (ICOS_prob - analytic_prob);
    }
  }
  printf("uniform diff squared: %f\n", sqrt(uni_diff_squared));
  printf("ICOS diff squared: %f\n", sqrt(ICOS_diff_squared));
}


// polyscope callback
void myCallback() {
    if (ImGui::BeginCombo("##combo1", all_polygons_current_item.c_str())){ // The second parameter is the label previewed before opening the combo.
      for (std::string tmp_str: all_polyhedra_items){ // This enables not having to have a const char* arr[]. Or maybe I'm just a noob.
          bool is_selected = (all_polygons_current_item == tmp_str.c_str()); // You can store your selection however you want, outside or inside your objects
          if (ImGui::Selectable(tmp_str.c_str(), is_selected)){ // selected smth
              polyscope::removeAllStructures();
              all_polygons_current_item = tmp_str;
              generate_polyhedron_example(all_polygons_current_item);
              auto t1 = clock();
              update_solver();
              printf("solver update time: %f\n", (clock() - t1)/(double)CLOCKS_PER_SEC);
              init_visuals();

              // visualize_colored_polyhedra();
              visualize_gauss_map();
              
              // G = find_center_of_mass(*forwardSolver->inputMesh, *forwardSolver->inputGeometry).first;
              // update_solver_and_boundaries();
              // boundary_builder->print_area_of_boundary_loops();
          }
          if (is_selected)
              ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
      }
      ImGui::EndCombo();
    }
    // if (ImGui::SliderInt("sim step count", &step_count, 1, 20));
    // if (ImGui::SliderFloat("orientation_vec X", &refresh_x, -1, 1) ||
    //     ImGui::SliderFloat("orientation_vec Y", &refresh_y, -1, 1) ||
    //     ImGui::SliderFloat("orientation_vec Z", &refresh_z, -1, 1)){
    //   refresh_orientation = Vector3({refresh_x, refresh_y, refresh_z}).normalize();
    //   old_g_vec = refresh_orientation;
    //   auto init_pc = polyscope::registerPointCloud("initial orientation", std::vector<Vector3>{refresh_orientation + vis_utils.center});
    //   init_pc->setPointColor({0.1,0.1,0.1});
    //   init_pc->setPointRadius(0.01, false);
    //   init_pc->setEnabled(true);
    // }
    // if (ImGui::Button("refresh")){
    //     old_g_vec = refresh_orientation;
    // }
    if (ImGui::InputInt("sample count", &total_samples));
    // if (ImGui::Button("build the raster image")){
    //   // draw the default polyhedra
    //   visualize_colored_polyhedra();
    //   visualize_gauss_map();
    //   build_raster_image();
    // }
    if (ImGui::Button("draw MS complex")){
      draw_stable_patches_on_gauss_map();
    }
    // if (ImGui::Checkbox("Save pos to file", &save_pos_to_file));
    if (ImGui::Checkbox("sampling ICOS", &ICOS_sampling));
    // if (ImGui::InputInt("max steps IPC", &max_steps_IPC));
    if (ImGui::Button("run IPC simulation")){
      if (ICOS_sampling)
        run_IPC_samples_ICOS("example.json", total_samples, max_steps_IPC);
      else
        run_IPC_samples_MCMC("example.json", total_samples, max_steps_IPC);
    }
}


int main(int argc, char* argv[])
{
#ifdef BT_USE_DOUBLE_PRECISION
	printf("BT_USE_DOUBLE_PRECISION\n");
#else
  printf("Single precision\n");
#endif
  polyscope::init();
  vis_utils = VisualUtils();
  generate_polyhedron_example(all_polygons_current_item);
  update_solver();
  init_visuals();
//   get_inetrial_rot("example.json");
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  // polyscope::view::upDir = polyscope::view::UpDir::YUp;
  // polyscope::options::groundPlaneHeightFactor = 1.; // adjust the plane height
  polyscope::state::userCallback = myCallback;


  // Give control to the polyscope gui
  polyscope::show();

  
  return EXIT_SUCCESS;	
}