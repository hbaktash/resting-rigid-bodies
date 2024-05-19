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

#include "LinearMath/btVector3.h"
#include "btBulletDynamicsCommon.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
// #include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
// #include "bullet3/examples/BasicExample.h"
#include "args/args.hxx"
#include "imgui.h"

#include "coloring.h"
#include "forward3D.h"
#include "mesh_factory.h"
#include "geometry_utils.h"
#include "bullet_sim.h"
#include "visual_utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// simulation stuff
PhysicsEnv* my_env;
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
std::vector<std::string> all_polyhedra_items = {std::string("tet"), std::string("tet2"), std::string("cube"), std::string("tilted cube"), std::string("sliced tet"), std::string("fox"), std::string("small_bunny"), std::string("bunnylp"), std::string("kitten"), std::string("double-torus"), std::string("knuckle_bone_real"),std::string("soccerball"), std::string("bunny"), std::string("gomboc"), std::string("dragon1"), std::string("dragon3"), std::string("mark_gomboc"), std::string("KnuckleboneDice"), std::string("Duende"), std::string("papa_noel"), std::string("reno"), std::string("baby_car"), std::string("rubberDuckie")};
std::string all_polygons_current_item = "rubberDuckie";
static const char* all_polygons_current_item_c_str = "rubberDuckie";

// GC stuff
std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr;
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;


// raster image stuff
FaceData<Vector3> face_colors;
int sample_count = 1e4;


// quasi static simulation stuff
Forward3DSolver* forwardSolver;
BoundaryBuilder *boundary_builder;
VisualUtils vis_utils;


//GM stuff
ManifoldSurfaceMesh* sphere_mesh;
VertexPositionGeometry* sphere_geometry;
bool gm_is_drawn = false;
bool draw_snail_trail = true;
bool save_pos_to_file = false;
Vector3 old_g_vec, new_g_vec;
int snail_trail_dummy_counter = 0;

void update_positions(){
    Vector<Vector3> new_positions = my_env->get_new_positions(forwardSolver->inputGeometry->inputVertexPositions.toVector());
    polyscope::getSurfaceMesh("my polyhedra")->updateVertexPositions(new_positions);
    polyscope::getSurfaceMesh("my polyhedra")->setTransparency(0.55);
    // // testing center of mass update
    // VertexData<Vector3> new_positions_vd(*forwardSolver->inputMesh);
    // new_positions_vd.fromVector(new_positions);
    // VertexPositionGeometry new_geom(*forwardSolver->inputMesh, vertex_data_to_matrix(new_positions_vd));
    // Vector3 new_center_of_mass = find_center_of_mass(*forwardSolver->inputMesh, new_geom).first;
    // std::cout << "new center of mass is: " << new_center_of_mass << "\n";
    // std::cout << "bullet current COM   : " << my_env->get_current_G() << "\n";
    // polyscope::registerPointCloud("new COM", std::vector<Vector3>{new_center_of_mass});
    // polyscope::registerPointCloud("env COM", std::vector<Vector3>{my_env->get_current_G()});
    // Forward3DSolver *tmp_solver = new Forward3DSolver(forwardSolver->inputMesh, &new_geom, my_env->get_current_G(), true);
    // tmp_solver->initialize_pre_computes();
    // BoundaryBuilder *tmp_builder = new BoundaryBuilder(tmp_solver);
    // tmp_builder->build_boundary_normals();
    // printf("Testing COM issue\n");
    // tmp_builder->print_area_of_boundary_loops();

    Vector<Vector3> new_hull_positions = my_env->get_new_positions(forwardSolver->hullGeometry->inputVertexPositions.toVector());
    VertexData<Vector3> new_hull_positions_vd(*my_env->mesh);
    new_hull_positions_vd.fromVector(new_hull_positions);
    Face touch_face =  my_env->get_touching_face(new_hull_positions_vd);
    printf("touching face after update is %d\n", touch_face.getIndex());

    if (save_pos_to_file){
      VertexPositionGeometry new_geom(*my_env->mesh, vertex_data_to_matrix(new_hull_positions_vd));
      writeSurfaceMesh(*my_env->mesh, new_geom, "../meshes/restPoses/"+all_polygons_current_item+"_rest.obj");
      printf(" current center of mass: \n");
      std::cout << my_env->get_current_G() << "\n";
    }
    // TODO: do this maybe instead of updating in the env _/_
    // psMesh->setTransform()
}


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
    std::unique_ptr<SurfaceMesh> nm_mesh_ptr;
    std::unique_ptr<VertexPositionGeometry> nm_geometry_ptr;
    std::string filename = "../meshes/"+poly_str+".obj";
    std::tie(nm_mesh_ptr, geometry_ptr) = readSurfaceMesh(filename);
    SurfaceMesh *nm_mesh = nm_mesh_ptr.release();
    mesh_ptr = nm_mesh->toManifoldMesh();
    mesh = mesh_ptr.release();
    geometry = geometry_ptr.release();
    preprocess_mesh(mesh, geometry, true);
  //   printf("baby car read\n");
  // }
  // else{
  //   std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(poly_str);
  //   mesh = mesh_ptr.release();
  //   geometry = geometry_ptr.release();
  //   preprocess_mesh(mesh, geometry, triangulate || std::strcmp(poly_str.c_str(), "gomboc") == 0, false, 1.);
  // }
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
  // first_time = true;
}

void initalize_env(std::string poly_str){
  // physics env
  my_env = new PhysicsEnv();
  my_env->init_physics();
  my_env->init_geometry(forwardSolver->hullMesh, forwardSolver->hullGeometry);
  my_env->add_ground(ground_box_y, ground_box_shape);
  my_env->add_object(G, Vector3({0,-1,0}));

  // polyscope
  initialize_vis(true);
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

void init_visuals(){
  // Register the mesh with polyscope
  polyscope::registerSurfaceMesh(
    "temp input mesh",
    geometry->inputVertexPositions + Vector3({1.4,2,1.4}), mesh->getFaceVertexList(),
    polyscopePermutations(*mesh));
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
  gm_is_drawn = true;
}

void draw_stable_patches_on_gauss_map(bool on_height_surface = false){
  if (!forwardSolver->updated)
    forwardSolver->initialize_pre_computes();
  // std::vector<Vector3> boundary_normals;
  auto net_pair = build_and_draw_stable_patches_on_gauss_map(boundary_builder,
                                                              vis_utils.center, vis_utils.gm_radi, vis_utils.arcs_seg_count, 
                                                              on_height_surface);
}

// // copying from main polyhedra exec; should proly move these to a sep module
// // TODO check the face normal flipped issue!!!
// void visualize_gauss_map(){
//   // draw the default polyhedra
//   visualize_colored_polyhedra();
//   // just draw the sphere next to the main surface
//   std::vector<Vector3> sphere_pos = {shift};
//   auto gauss_map_pc = polyscope::registerPointCloud("Gauss Map", sphere_pos);
//   gauss_map_pc->setPointColor({0.74,0.7,0.9});
//   gauss_map_pc->setPointRadius(gm_radi, false);
//   gauss_map_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

//   // point cloud for face normals; normals should be outwards
//   std::vector<Vector3> face_normal_points, face_normal_colors;
//   for (Face f: forwardSolver->hullMesh->faces()){
//     Vector3 normal_pos_on_gm = forwardSolver->hullGeometry->faceNormal(f) + shift;
//     face_normal_points.push_back(normal_pos_on_gm);
//     face_normal_colors.push_back(face_colors[f]);
//   }
//   auto face_normals_pc = polyscope::registerPointCloud("Face Normals", face_normal_points);
//   face_normals_pc->setPointRadius(face_normal_vertex_gm_radi, false);
//   face_normals_pc->setPointColor({0.9,0.9,0.9});
//   face_normals_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);
//   polyscope::PointCloudColorQuantity* face_normals_pc_colors = face_normals_pc->addColorQuantity("random color", face_normal_colors);
//   face_normals_pc_colors->setEnabled(true);
  
//   // arcs for edge-normals set
//   std::vector<std::vector<size_t>> dummy_face{{0,0,0}};
//   //    dummy mesh to add curves to
//   polyscope::registerSurfaceMesh(
//       "dummy mesh for gauss map arcs",
//       geometry->inputVertexPositions, dummy_face);
//   //    add arc per edge
//   for (Edge e: forwardSolver->hullMesh->edges()){
//     Face f1 = e.halfedge().face(),
//          f2 = e.halfedge().twin().face();
//     Vector3 n1 = forwardSolver->hullGeometry->faceNormal(f1),
//             n2 = forwardSolver->hullGeometry->faceNormal(f2);
//     // draw with polyscope
//     draw_arc_on_sphere(n1, n2, shift, gm_radi, arcs_seg_count, e.getIndex());
//   }
// }


// sample and raster
void build_raster_image(){
  FaceData<std::vector<Vector3>> face_samples(*my_env->mesh);
  std::vector<Vector3> final_orientations;
  std::vector<Vector3> falsh;
  int total_invalids = 0, total_samples = 0;
  auto t1 = clock();
  double avg_time = 0.;
  for (int i = 0; i < sample_count; i++){
    Vector3 random_orientation = {randomReal(-1,1), randomReal(-1,1), randomReal(-1,1)};
    if (random_orientation.norm() <= 1){
      auto t2 = clock();
      total_samples++;
      if (i % 2000 == 0)
        printf("$$$ at sample %d\n", i);
      random_orientation /= norm(random_orientation);
      my_env->refresh(G, random_orientation);
      Face touching_face = my_env->final_stable_face(); // Invalid if not close to a face normal; wtf?
      avg_time += (double)((clock() - t2)/(double)CLOCKS_PER_SEC);
      final_orientations.push_back(my_env->get_current_orientation() + vis_utils.center);
      if (touching_face.getIndex() == INVALID_IND){
        total_invalids++;
        falsh.push_back(random_orientation + vis_utils.center); //  shift for visualization
        continue;
      }
      face_samples[touching_face].push_back(random_orientation);
    }
  }
  printf(" --- sampling time: %f\n", (double)(clock() - t1));
  printf(" ### total invalid faces: %d/%d\n", total_invalids, total_samples);
  std::vector<Vector3> raster_positions,
                       raster_colors;
  
  forwardSolver->build_face_last_faces();
  FaceData<double> accum_face_areas(*forwardSolver->hullMesh, 0.);
  printf("empirical probs:\n");
  size_t unstable_faces = 0;
  for (Face f: my_env->mesh->faces()){
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
  for (Face f: my_env->mesh->faces()){
    if (forwardSolver->face_is_stable(f)){
      printf("  f %d -> emp: %f    ----   MS: %f \n", f.getIndex(), accum_face_areas[f], boundary_builder->face_region_area[f]/(4.*PI));
    }
  }
  //
  printf(" $$$ KL divergence! $$$\n");
  double kl_divergence = 0.;
  for (Face f: my_env->mesh->faces()){
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


// void test_static_dice_pipeline(){
//   // forwardSolver->set_uniform_G();
//   ManifoldSurfaceMesh *hull_mesh = new ManifoldSurfaceMesh(forwardSolver->hullMesh->getFaceVertexList());
//   VertexPositionGeometry *hull_geo = new VertexPositionGeometry(*hull_mesh, vertex_data_to_matrix(forwardSolver->hullGeometry->inputVertexPositions));
//   Forward3DSolver *tmp_solver = new Forward3DSolver(hull_mesh, hull_geo, G, false);
  
//   auto GV_pair = find_center_of_mass(*forwardSolver->inputMesh, *forwardSolver->inputGeometry);
//   tmp_solver->set_G(GV_pair.first);
//   tmp_solver->volume = GV_pair.second;
  
//   printf("initalize precomputes\n");
//   tmp_solver->initialize_pre_computes();

//   BoundaryBuilder *tmp_bnd_builder = new BoundaryBuilder(tmp_solver);
//   tmp_bnd_builder->build_boundary_normals(); // (poses_ad, G_ad, ) // autodiff; generate_gradients = true
  
//   Vector3 GG = tmp_solver->get_G();
//   tmp_solver->build_face_last_faces();
//   printf("testing static dice E\n");
//   tmp_bnd_builder->dice_energy(vertex_data_to_matrix(tmp_solver->hullGeometry->inputVertexPositions),
//                                Eigen::Vector3d({GG.x, GG.y, GG.z}), *hull_mesh,
//                                tmp_bnd_builder->find_terminal_edges(), tmp_solver->face_last_face, tmp_solver->vertex_is_stabilizable, tmp_solver->edge_next_vertex,
//                                6);
  
//   // update_visuals_with_G(tmp_solver, tmp_bnd_builder);
// }


void find_best_duck_G(){
  double initial_low_face_probs = 0.;
  for (Face f: forwardSolver->hullMesh->faces()){
    if ( (forwardSolver->hullGeometry->faceNormal(f) - Vector3({0,-1.,0})).norm() < 0.2){
      initial_low_face_probs += boundary_builder->face_region_area[f]/(4.*PI);
    }
  }
  int resolution_x = 10,
      resolution_z = 40;
  double h_y = 0.;
  double highest_stand_prob = 0.0;
  Vector3 best_G;
  for (int i = 0; i < resolution_x; i++){
    double h_x = -0.2 + 0.4*(double)i/(double)resolution_x;
    for (int j = 0; j < resolution_z; j++){
      double h_z = -0.5 + 1.*(double)j/(double)resolution_z;
      Vector3 tmp_G({h_x, h_y, h_z});
      forwardSolver->set_G(tmp_G);
      forwardSolver->initialize_pre_computes();
      boundary_builder = new BoundaryBuilder(forwardSolver);
      boundary_builder->build_boundary_normals();
      double low_face_probs = 0.;
      for (Face f: forwardSolver->hullMesh->faces()){
        if ( (forwardSolver->hullGeometry->faceNormal(f) - Vector3({0,-1.,0})).norm() < 0.2){
          low_face_probs += boundary_builder->face_region_area[f]/(4.*PI);
        }
      }
      if (low_face_probs > highest_stand_prob){
        highest_stand_prob = low_face_probs;
        best_G = tmp_G;
        printf("current best prob  : %f\n", highest_stand_prob);
      }
    }
  }
  forwardSolver->set_G(best_G);
  forwardSolver->initialize_pre_computes();
  boundary_builder = new BoundaryBuilder(forwardSolver);
  boundary_builder->build_boundary_normals();
  
  printf("best G is: %f, %f, %f\n", best_G.x, best_G.y, best_G.z);
  printf("with best prob        : %f\n", highest_stand_prob);   
  printf("initial low face probs: %f\n", initial_low_face_probs);
  draw_stable_patches_on_gauss_map();
  auto G_pc = polyscope::registerPointCloud("best COM", std::vector<Vector3>{best_G});
  G_pc->setPointColor({0.1,0.1,0.1});
  G_pc->setPointRadius(0.03, false);
  G_pc->setEnabled(true);
  visualize_colored_polyhedra();
}


void find_best_bunny_G(){
  size_t side_count = 6;
  double old_dice_E = boundary_builder->get_fair_dice_energy(side_count);
  int resolution_x = 20,
      resolution_z = 20,
      resolution_y = 20;
  // double x0 = 0.23,
  //        y0 = 0.23,
  //        z0 = 0.23;
  double x0 = 0.022,
         y0 = 0.14,
         z0 = 0.2;
  double lowest_dice_E = old_dice_E;
  Vector3 best_G;
  for (int i = 0; i < resolution_x; i++){
    printf(" at i: %d\n", i);
    double h_x = -x0 + x0*0.5*(double)i/(double)resolution_x;
    for (int k = 0; k < resolution_y; k++){
      double h_y = -y0 + y0*0.5*(double)k/(double)resolution_y;
      for (int j = 0; j < resolution_z; j++){
        double h_z = -z0 + z0*0.5*(double)j/(double)resolution_z;
        Vector3 tmp_G({h_x, h_y, h_z});
        forwardSolver->set_G(tmp_G);
        forwardSolver->initialize_pre_computes();
        boundary_builder = new BoundaryBuilder(forwardSolver);
        boundary_builder->build_boundary_normals();
        double dice_E = boundary_builder->get_fair_dice_energy(side_count);
        if (dice_E < lowest_dice_E){
          lowest_dice_E = dice_E;
          best_G = tmp_G;
          printf("current dice E  : %f\n", lowest_dice_E);
        }
      }
    }
  }
  forwardSolver->set_G(best_G);
  forwardSolver->initialize_pre_computes();
  boundary_builder = new BoundaryBuilder(forwardSolver);
  boundary_builder->build_boundary_normals();
  boundary_builder->print_area_of_boundary_loops();
  printf("best G is: %f, %f, %f\n", best_G.x, best_G.y, best_G.z);
  printf("with best Dice E        : %f\n", lowest_dice_E);   
  printf("initial dice E: %f\n", old_dice_E);
  draw_stable_patches_on_gauss_map();
  auto G_pc = polyscope::registerPointCloud("best COM", std::vector<Vector3>{best_G});
  G_pc->setPointColor({0.1,0.9,0.1});
  G_pc->setPointRadius(0.2, false);
  G_pc->setEnabled(true);
  visualize_colored_polyhedra();
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
              initalize_env(all_polygons_current_item);
              init_visuals();

              visualize_colored_polyhedra();
              visualize_gauss_map();
              // //
              // G = find_center_of_mass(*forwardSolver->inputMesh, *forwardSolver->inputGeometry).first;
              // update_solver_and_boundaries();
              // boundary_builder->print_area_of_boundary_loops();
          }
          if (is_selected)
              ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
      }
      ImGui::EndCombo();
  }
    if (ImGui::Button("take simulation step")){
        my_env->take_step(step_count, step_size);
        if(draw_snail_trail){
          new_g_vec = my_env->get_current_orientation();
          if (norm(old_g_vec-new_g_vec) != 0.){ // proly dont have to use tol
            draw_arc_on_sphere(old_g_vec, new_g_vec, vis_utils.center,
                               vis_utils.gm_radi, vis_utils.arcs_seg_count, 200 , //+snail_trail_dummy_counter, 
                               1.2, {0.9,0.8,0.1});
            snail_trail_dummy_counter++;
          }
          old_g_vec = new_g_vec;
        }
        update_positions();
    }
    if(ImGui::SliderFloat("sim step size", &step_size, 0.0001, 0.1)) my_env->default_step_size = step_size;
    if(ImGui::SliderInt("sim step count", &step_count, 1, 20));
    
    if (ImGui::Button("fast forward to stable state")){
      auto t1 = clock();
      Face touching_face = my_env->final_stable_face(draw_snail_trail);
      printf("time to stable state: %f\n", (clock() - t1)/(double)CLOCKS_PER_SEC);
      if (draw_snail_trail){
        std::vector<Vector3> snail_trail = my_env->orientation_trail;
        // std::vector<Vector3> snail_trail;
        // for (size_t i = 0; i < orig_snail_trail.size(); i++){
        //   if (i % 2 == 0)
        //     snail_trail.push_back(orig_snail_trail[i]);
        // }
        for (Vector3 &v: snail_trail)
          v += vis_utils.center;
        auto trail_pc = polyscope::registerPointCloud("bullet pc trail", snail_trail);
        glm::vec3 init_color = {0.8,0.8,0.2};
        std::vector<glm::vec3> colors;
        for (size_t i = 0; i < snail_trail.size(); i++){
          glm::vec3 tmp_color = init_color;
          tmp_color = tmp_color * (1.f - (i/(float)snail_trail.size()));
          colors.push_back(tmp_color);
        }
        trail_pc->addColorQuantity("time", colors)->setEnabled(true);
        trail_pc->setEnabled(true);

        // draw_trail_on_gm(snail_trail, {0.9,0.8,0.1}, "bullet trail",1.0, true);
      }
      printf("final touching face is %d\n", touching_face.getIndex());
      if (touching_face.getIndex() != INVALID_IND)
        std::cout << "face normal is "<< forwardSolver->hullGeometry->faceNormal(touching_face)<< "\n";
      update_positions();
    }
    if (ImGui::Checkbox("snail trail", &draw_snail_trail));
    if (ImGui::Button("Build quasistatic snail trail")){
      auto t1 = clock();
      std::vector<Vector3> snail_trail = forwardSolver->snail_trail_log(refresh_orientation);
      printf(" Quasistatic trail time: %f\n", (clock() - t1)/(double)CLOCKS_PER_SEC);
      draw_trail_on_gm(snail_trail, {0.7,0.1,0.8}, "quasi-static trail",1.);
    }
    if (ImGui::SliderFloat("orientation_vec X", &refresh_x, -1, 1) ||
        ImGui::SliderFloat("orientation_vec Y", &refresh_y, -1, 1) ||
        ImGui::SliderFloat("orientation_vec Z", &refresh_z, -1, 1)){
      refresh_orientation = Vector3({refresh_x, refresh_y, refresh_z}).normalize();
      old_g_vec = refresh_orientation;
      auto init_pc = polyscope::registerPointCloud("initial orientation", std::vector<Vector3>{refresh_orientation + vis_utils.center});
      init_pc->setPointColor({0.1,0.1,0.1});
      init_pc->setPointRadius(0.01, false);
      init_pc->setEnabled(true);
    }
    if (ImGui::Button("refresh")){
        my_env->refresh(G, refresh_orientation);
        old_g_vec = refresh_orientation;
    }
    if (ImGui::InputInt("sample count", &sample_count));
    if (ImGui::Button("build the raster image")){
      
      // draw the default polyhedra
      visualize_colored_polyhedra();
      visualize_gauss_map();
      build_raster_image();
    }
    if (ImGui::Button("draw MS complex")){
      draw_stable_patches_on_gauss_map();
    }
    if (ImGui::Checkbox("Save pos to file", &save_pos_to_file));

    if (ImGui::SliderFloat("G Y", &G_y, -0.5, 0.5) ||
        ImGui::SliderFloat("G X", &G_x, -0.5, 0.5) ||
        ImGui::SliderFloat("G Z", &G_z, -0.5, 0.5)){
      G = Vector3({G_x, G_y, G_z});
      forwardSolver->set_G(G);
      forwardSolver->initialize_pre_computes();
      boundary_builder = new BoundaryBuilder(forwardSolver);
      boundary_builder->build_boundary_normals();
      draw_stable_patches_on_gauss_map();
      auto G_pc = polyscope::registerPointCloud("COM", std::vector<Vector3>{G});
      G_pc->setPointColor({0.1,0.1,0.1});
      G_pc->setPointRadius(0.03, false);
      G_pc->setEnabled(true);
      visualize_colored_polyhedra();
      if (forwardSolver->hullMesh->nFaces() > 300){ // for baby car
        double p1 = boundary_builder->face_region_area[forwardSolver->hullMesh->face(298)]/(4.*PI);
        double p2 = boundary_builder->face_region_area[forwardSolver->hullMesh->face(297)]/(4.*PI);
        if (p1 > 0 && p2 == 0)
          printf(" face 298 prob: %f\n", p1);
        else if (p2 > 0 && p1 == 0)
          printf(" face 297 prob: %f\n", p2);
        else
          printf(" face 297+298 prob: %f\n", p1+p2);
      }
    }
    if (ImGui::Button("highlight low faces")){
      face_colors = FaceData<Vector3>(*forwardSolver->hullMesh, default_face_color);
      for (Face f: forwardSolver->hullMesh->faces()){
        if ( (forwardSolver->hullGeometry->faceNormal(f) - Vector3({0,-1.,0})).norm() < 0.2){
          face_colors[f] = {0.,0.,0.};
        }
      }
      polyscope::SurfaceFaceColorQuantity *faceQnty2 = polyscope::getSurfaceMesh("init hull mesh")->addFaceColorQuantity("random face colors2", face_colors);
      faceQnty2->setEnabled(true);
    }
    if (ImGui::Button("find best duck G")){
      find_best_duck_G();
    }
    if (ImGui::Button("find best Bunny G")){
      find_best_bunny_G();
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

  // physics env
  initalize_env(all_polygons_current_item);
  
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  // polyscope::view::upDir = polyscope::view::UpDir::YUp;
  // polyscope::options::groundPlaneHeightFactor = 1.; // adjust the plane height
  polyscope::state::userCallback = myCallback;


  // Give control to the polyscope gui
  polyscope::show();

  
  return EXIT_SUCCESS;
	
}
