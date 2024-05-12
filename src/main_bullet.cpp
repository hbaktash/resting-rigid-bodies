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
float step_size = 0.006, // 0.016
      refresh_x, refresh_y, refresh_z;
int step_count = 1;
Vector3 G;
Vector3 refresh_orientation({0,-1,0});

// stuff for Gauss map
float face_normal_vertex_gm_radi = 0.03,
      gm_distance = 2.1,
      gm_radi = 1.;
int arcs_seg_count = 13;
Vector3 shift = {0., gm_distance , 0.},
        colored_shift = {gm_distance, gm_distance , 0.};
float arc_curve_radi = 0.01;


double ground_box_y = -2;
Vector3 ground_box_shape({10,1,10});

Vector3 default_face_color({0.99,0.99,0.99});

// example choice
std::vector<std::string> all_polyhedra_items = {std::string("tet"), std::string("tet2"), std::string("cube"), std::string("tilted cube"), std::string("sliced tet"), std::string("Conway spiral 4"), std::string("oloid"), std::string("fox"), std::string("small_bunny"), std::string("bunnylp"), std::string("kitten"), std::string("double-torus"), std::string("soccerball"), std::string("cowhead"), std::string("bunny"), std::string("gomboc"), std::string("mark_gomboc")};
std::string all_polygons_current_item = "tet";
static const char* all_polygons_current_item_c_str = "tet";

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


void update_positions(){
    Vector<Vector3> new_positions = my_env->get_new_positions(forwardSolver->inputGeometry->inputVertexPositions.toVector());
    polyscope::getSurfaceMesh("my polyhedra")->updateVertexPositions(new_positions);
    
    Vector<Vector3> new_hull_positions = my_env->get_new_positions(forwardSolver->hullGeometry->inputVertexPositions.toVector());
    VertexData<Vector3> new_hull_positions_vd(*my_env->mesh);
    new_hull_positions_vd.fromVector(new_hull_positions);
    Face touch_face =  my_env->get_touching_face(new_hull_positions_vd);
    printf("touching face after update is %d\n", touch_face.getIndex());
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
  std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(poly_str);
  mesh = mesh_ptr.release();
  geometry = geometry_ptr.release();
  preprocess_mesh(mesh, geometry, triangulate || std::strcmp(poly_str.c_str(), "gomboc") == 0, false, 1.);
  G = find_center_of_mass(*mesh, *geometry).first;
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
}


void update_solver(){
  //assuming convex input here
  forwardSolver = new Forward3DSolver(mesh, geometry, G, true);
  forwardSolver->initialize_pre_computes();
  initialize_boundary_builder();
  vis_utils.forwardSolver = forwardSolver;
}

void init_visuals(){
  // Register the mesh with polyscope
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
  face_colors = generate_random_colors(fwd_solver->hullMesh);
  
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
  int total_invalids = 0, total_samples = 0;
  for (int i = 0; i < sample_count; i++){
    Vector3 random_orientation = {randomReal(-1,1), randomReal(-1,1), randomReal(-1,1)};
    if (random_orientation.norm() <= 1){
      total_samples++;
      if (i % 2000 == 0)
        printf("$$$ at sample %d\n", i);
      random_orientation /= norm(random_orientation);
      my_env->refresh(G, random_orientation);
      Face touching_face = my_env->final_stable_face();
      if (touching_face.getIndex() == INVALID_IND){
        total_invalids++;
        continue;
      }
      // if (touching_face.getIndex() == 0){
        // std::cout<<" -- face 0 vec: "<< random_orientation << "\n";
      // }
      // printf("touching face is %d\n", touching_face.getIndex());
      face_samples[touching_face].push_back(random_orientation);
    }
  }
  printf(" ### total invalid faces: %d/%d\n", total_invalids, total_samples);
  std::vector<Vector3> raster_positions,
                       raster_colors;
  
  forwardSolver->build_face_last_faces();
  FaceData<double> accum_facee_areas(*forwardSolver->hullMesh, 0.);
  printf("empirical probs:\n");
  for (Face f: my_env->mesh->faces()){
    std::vector<Vector3> tmp_points = face_samples[f];
    double prob = tmp_points.size()/(double)(total_samples - total_invalids);
    accum_facee_areas[forwardSolver->face_last_face[f]] += prob;
    if(tmp_points.size() != 0) printf(" --- f %d: %f\n", f.getIndex(), prob);
    for (Vector3 tmp_p: tmp_points){
      raster_positions.push_back(tmp_p + shift);
      raster_colors.push_back(face_colors[forwardSolver->face_last_face[f]]);
      // std::cout<< "tmp color is: " << face_colors[f] << "\n";
    }
  }
  printf("accumulated probs:\n");
  for (Face f: my_env->mesh->faces()){
    if (forwardSolver->face_is_stable(f)){
      printf("  f %d: %f\n", f.getIndex(), accum_facee_areas[f]);
    }
  }
  printf("morse complex probs:\n");
  for (Face f: my_env->mesh->faces()){
    if (forwardSolver->face_is_stable(f)){
      printf("  f %d: %f\n", f.getIndex(), boundary_builder->face_region_area[f]/(4.*PI));
    }
  }
  // boundary_builder->print_area_of_boundary_loops();
  auto raster_pc = polyscope::registerPointCloud("raster point cloud", raster_positions);
  polyscope::PointCloudColorQuantity* pc_col_quant = raster_pc->addColorQuantity("random color", raster_colors);
  pc_col_quant->setEnabled(true);
}


// polyscope callback
void myCallback() {
    if (ImGui::BeginCombo("##combo1", all_polygons_current_item.c_str())){ // The second parameter is the label previewed before opening the combo.
      for (std::string tmp_str: all_polyhedra_items){ // This enables not having to have a const char* arr[]. Or maybe I'm just a noob.
          bool is_selected = (all_polygons_current_item == tmp_str.c_str()); // You can store your selection however you want, outside or inside your objects
          if (ImGui::Selectable(tmp_str.c_str(), is_selected)){ // selected smth
              all_polygons_current_item = tmp_str;
              generate_polyhedron_example(all_polygons_current_item);
              update_solver();
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
        update_positions();
    }
    if(ImGui::SliderFloat("sim step size", &step_size, 0.005, 0.1)) my_env->default_step_size = step_size;
    if(ImGui::SliderInt("sim step count", &step_count, 1, 20));
    
    if (ImGui::Button("fast forward to stable state")){
        Face touching_face = my_env->final_stable_face();
        printf("final touching face is %d\n", touching_face.getIndex());
        if (touching_face.getIndex() != INVALID_IND)
          std::cout << "face normal is "<< forwardSolver->hullGeometry->faceNormal(touching_face)<< "\n";
        update_positions();
    }
    if (ImGui::InputFloat("orientation_vec X", &refresh_x) ||
        ImGui::InputFloat("orientation_vec Y", &refresh_y) ||
        ImGui::InputFloat("orientation_vec Z", &refresh_z)){
      refresh_orientation = {refresh_x, refresh_y, refresh_z};
    }
    if (ImGui::Button("refresh")){
        my_env->refresh(G, refresh_orientation);
    }
    if (ImGui::Button("Show Gauss Map") || 
        ImGui::SliderInt("seg count for arcs", &arcs_seg_count, 1, 100)||
        ImGui::SliderFloat("arc curve radi", &arc_curve_radi, 0., 0.04)||
        ImGui::SliderFloat("face normal vertex radi", &face_normal_vertex_gm_radi, 0., 0.04)){///face_normal_vertex_gm_radi
            
          // draw the default polyhedra
          visualize_colored_polyhedra();
          visualize_gauss_map();
    }
    
    if (ImGui::Button("draw samples on Gauss map")){
        Face touching_face = my_env->final_stable_face();
        printf("final touching face is %d\n", touching_face.getIndex());
        std::cout << "face normal is "<< forwardSolver->hullGeometry->faceNormal(touching_face)<< "\n";
        update_positions();
    }
    if (ImGui::InputInt("sample count", &sample_count));
    if (ImGui::Button("build the raster image")){
      
      // draw the default polyhedra
      visualize_colored_polyhedra();
      visualize_gauss_map();
      build_raster_image();
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
  vis_utils.gm_distance = gm_distance;
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
