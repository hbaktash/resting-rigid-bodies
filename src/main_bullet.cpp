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
#include "mesh_factory.h"
#include "geometry_utils.h"
#include "bullet_sim.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// simulation stuff
PhysicsEnv* my_env;
float step_size = 0.016,
      refresh_x, refresh_y, refresh_z;
int step_count = 1;
Vector3 G;
Vector3 refresh_g_vec({0,-1,0});

// stuff for Gauss map
float arc_curve_radi = 0.01,
      face_normal_vertex_gm_radi = 0.03,
      gm_distance = 2.1,
      gm_radi = 1.;
int arcs_seg_count = 13;
Vector3 shift = {0., 0. , gm_distance},
        colored_shift = {0., 0. , -gm_distance};

double ground_box_y = -3;
Vector3 ground_box_shape({10,1,10});

polyscope::SlicePlane* psPlane;



// example choice
std::vector<std::string> all_polyhedra_items = {std::string("tet"), std::string("tet2"), std::string("sliced tet")};
std::string all_polygons_current_item = "tet";
static const char* all_polygons_current_item_c_str = "tet";

// GC stuff
std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr;
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

// polyscope stuff
polyscope::SurfaceMesh *psMesh, *coloredPsMesh, *dummy_psMesh2;
polyscope::PointCloud *gauss_map_pc, *face_normals_pc, *raster_pc;

// raster image stuff
FaceData<Vector3> face_colors;
int sample_count = 5e4;


void update_positions(){
    VertexData<Vector3> new_positions = my_env->get_new_positions();
    psMesh->updateVertexPositions(new_positions.raw());
    // TODO: do this maybe instead of updating in the env _/_
    // psMesh->setTransform()
}


void initialize_vis(bool with_plane = true){
    psMesh = polyscope::registerSurfaceMesh("my polyhedra", geometry->inputVertexPositions, mesh->getFaceVertexList());
    // ground plane on Polyscope has a weird height setting (scaled factor..)
    if (with_plane){
      psPlane = polyscope::addSceneSlicePlane();
      psPlane->setDrawPlane(true);  // render the semi-transparent gridded plane
      psPlane->setDrawWidget(false);
      psPlane->setPose(glm::vec3{0., ground_box_y + 1, 0.}, glm::vec3{0., 1., 0.});
    }
}


void initalize_mesh_and_env(std::string poly_str){
  // mesh stuff
  std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(poly_str);
  mesh = mesh_ptr.release(); 
  geometry = geometry_ptr.release();
  center_and_normalize(mesh, geometry);
  G = find_center_of_mass(*mesh, *geometry);

  // physics env stuff
  my_env = new PhysicsEnv();
  my_env->init_physics();
  my_env->init_geometry(mesh, geometry);
  my_env->add_ground(ground_box_y, ground_box_shape);
  my_env->add_object(G, Vector3({0,-1,0}));

  // polyscope
  initialize_vis(false);
}


// another polyhedra for the sake of a good colored raster image
void visualize_colored_polyhedra(){
  VertexData<Vector3> shifted_positions(*mesh);
  for (Vertex v: mesh->vertices()){
    shifted_positions[v] = geometry->inputVertexPositions[v] + colored_shift;
  }
  coloredPsMesh = polyscope::registerSurfaceMesh("colored polyhedra", shifted_positions, mesh->getFaceVertexList());

  // generate random colors and color the faces
  face_colors = generate_random_colors(mesh);
  polyscope::SurfaceFaceColorQuantity *faceQnty = coloredPsMesh->addFaceColorQuantity("random face colors", face_colors);
  faceQnty->setEnabled(true);

  // add colors to the original polyhedra as well
  polyscope::SurfaceFaceColorQuantity *faceQnty2 = psMesh->addFaceColorQuantity("random face colors2", face_colors);
  faceQnty2->setEnabled(true);
}

// draw an arc connecting two points on the sphere; for Gauss map purposes
void draw_arc_on_sphere(Vector3 p1, Vector3 p2, Vector3 center, double radius, size_t seg_count, size_t edge_ind){
  // p1, p2 just represent normal vectors
  if (norm(p1) > 1.01 || norm(p2) > 1.01)
    polyscope::warning("wtf? p1, p2 norm larger than 1");
  
  std::vector<std::array<size_t, 2>> edgeInds;
  std::vector<Vector3> positions;
  double sqrt_radi = sqrt(radius);
  // walk on p1-p2 segment
  Vector3 curr_point = p1,
          forward_vec = (p2-p1)/(double)seg_count;
  Vector3 next_point = curr_point + forward_vec;
  Vector3 curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center ,
          next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
  positions.push_back(curr_point_on_sphere);
  for (size_t i = 0; i < seg_count; i++){
    // add to positions list
    curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center ,
    next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
    positions.push_back(next_point_on_sphere);
    // add segment indices
    edgeInds.push_back({i, i+1});

    // update points
    curr_point = next_point;
    next_point += forward_vec;
  }
  polyscope::SurfaceGraphQuantity* psArcCurve = dummy_psMesh2->addSurfaceGraphQuantity("Arc curve " + std::to_string(edge_ind), positions, edgeInds);
  psArcCurve->setRadius(arc_curve_radi, false);
  if (edge_ind < 100)
    psArcCurve->setColor({0.03, 0.03, 0.03});
  else
    psArcCurve->setColor({0.05, 0.5, 0.5});
  psArcCurve->setEnabled(true);
}

// copying from main polyhedra exec; should proly move these to a sep module
// TODO check the face normal flipped issue!!!
void visualize_gauss_map(){
  // draw the default polyhedra
  visualize_colored_polyhedra();
  // just draw the sphere next to the main surface
  std::vector<Vector3> sphere_pos = {shift};
  gauss_map_pc = polyscope::registerPointCloud("Gauss Map", sphere_pos);
  gauss_map_pc->setPointColor({0.74,0.7,0.9});
  gauss_map_pc->setPointRadius(gm_radi, false);
  gauss_map_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

  // point cloud for face normals; normals should be outwards
  std::vector<Vector3> face_normal_points, face_normal_colors;
  for (Face f: mesh->faces()){
    Vector3 normal_pos_on_gm = geometry->faceNormal(f) + shift;
    face_normal_points.push_back(normal_pos_on_gm);
    face_normal_colors.push_back(face_colors[f]);
  }
  face_normals_pc = polyscope::registerPointCloud("Face Normals", face_normal_points);
  face_normals_pc->setPointRadius(face_normal_vertex_gm_radi, false);
  face_normals_pc->setPointColor({0.9,0.9,0.9});
  face_normals_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  polyscope::PointCloudColorQuantity* face_normals_pc_colors = face_normals_pc->addColorQuantity("random color", face_normal_colors);
  face_normals_pc_colors->setEnabled(true);
  
  // arcs for edge-normals set
  std::vector<std::vector<size_t>> dummy_face{{0,0,0}};
  //    dummy mesh to add curves to
  dummy_psMesh2 = polyscope::registerSurfaceMesh(
      "dummy mesh for gauss map arcs",
      geometry->inputVertexPositions, dummy_face);
  //    add arc per edge
  for (Edge e: mesh->edges()){
    Face f1 = e.halfedge().face(),
         f2 = e.halfedge().twin().face();
    Vector3 n1 = geometry->faceNormal(f1),
            n2 = geometry->faceNormal(f2);
    // draw with polyscope
    draw_arc_on_sphere(n1, n2, shift, gm_radi, arcs_seg_count, e.getIndex());
  }
}


// sample and raster
void build_raster_image(){
  FaceData<std::vector<Vector3>> face_samples(*mesh);
  int total_invalids = 0;
  for (int i = 0; i < sample_count; i++){
    Vector3 random_g_vec = {randomReal(-1,1), randomReal(-1,1), randomReal(-1,1)};
    if (random_g_vec.norm() <= 1){
      if (i % 2000 == 0)
        printf("$$$ at sample %d\n", i);
      random_g_vec /= norm(random_g_vec);
      my_env->refresh(G, random_g_vec);
      Face touching_face = my_env->final_stable_face(Vector3({0,-1,0}));
      if (touching_face.getIndex() == INVALID_IND){
        total_invalids++;
        continue;
      }
      // if (touching_face.getIndex() == 0){
        // std::cout<<" -- face 0 vec: "<< random_g_vec<< "\n";
      // }
      face_samples[touching_face].push_back(random_g_vec);
    }
  }
  printf(" ### total invalid faces: %d\n", total_invalids);
  std::vector<Vector3> raster_positions,
                       raster_colors;
  for (Face f: mesh->faces()){
    std::vector<Vector3> tmp_points = face_samples[f];
    for (Vector3 tmp_p: tmp_points){
      raster_positions.push_back(tmp_p + shift);
      raster_colors.push_back(face_colors[f]);
      // std::cout<< "tmp color is: " << face_colors[f] << "\n";
    }
  }
  raster_pc = polyscope::registerPointCloud("raster point cloud", raster_positions);
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
              initalize_mesh_and_env(all_polygons_current_item);
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
    if(ImGui::SliderFloat("sim step size", &step_size, 0.01, 0.1));
    if(ImGui::SliderInt("sim step count", &step_count, 1, 20));
    
    if (ImGui::Button("fast forward to stable state")){
        Face touching_face = my_env->final_stable_face(Vector3({0,-1,0}));
        printf("final touching face is %d\n", touching_face.getIndex());
        std::cout << "face normal is "<< geometry->faceNormal(touching_face)<< "\n";
        update_positions();
    }
    if (ImGui::InputFloat("g_vec X", &refresh_x) ||
        ImGui::InputFloat("g_vec Y", &refresh_y) ||
        ImGui::InputFloat("g_vec Z", &refresh_z)){
      refresh_g_vec = {refresh_x, refresh_y, refresh_z};
    }
    if (ImGui::Button("refresh")){
        my_env->refresh(G, refresh_g_vec);
    }
    
    if (ImGui::Button("Show Gauss Map") || 
        ImGui::SliderInt("seg count for arcs", &arcs_seg_count, 1, 100)||
        ImGui::SliderFloat("arc curve radi", &arc_curve_radi, 0., 0.04)||
        ImGui::SliderFloat("face normal vertex radi", &face_normal_vertex_gm_radi, 0., 0.04)){///face_normal_vertex_gm_radi
            visualize_gauss_map();
    }
    
    if (ImGui::Button("draw samples on Gauss map")){
        Face touching_face = my_env->final_stable_face(Vector3({0,-1,0}));
        printf("final touching face is %d\n", touching_face.getIndex());
        std::cout << "face normal is "<< geometry->faceNormal(touching_face)<< "\n";
        update_positions();
    }
    
    if (ImGui::Button("build the raster image")){
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
  std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra("tet");
  mesh = mesh_ptr.release(); 
  geometry = geometry_ptr.release();
  center_and_normalize(mesh, geometry);
  G = find_center_of_mass(*mesh, *geometry);

  my_env = new PhysicsEnv();
  my_env->init_physics();
  my_env->init_geometry(mesh, geometry);

  ///-----initialization_end-----
  my_env->add_ground(ground_box_y, ground_box_shape);
  my_env->add_object(G, Vector3({0,-1,0}));

  polyscope::init();
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  // polyscope::view::upDir = polyscope::view::UpDir::YUp;
  // polyscope::options::groundPlaneHeightFactor = 1.; // adjust the plane height
  polyscope::state::userCallback = myCallback;
  initialize_vis();


  // Give control to the polyscope gui
  polyscope::show();

  
  return EXIT_SUCCESS;
	
}
