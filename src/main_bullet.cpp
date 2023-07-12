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

#include "mesh_factory.h"
#include "geometry_utils.h"
#include "bullet_sim.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// simulation stuff
PhysicsEnv* my_env;
float step_size = 0.03;
int step_count = 1;

// stuff for Gauss map
float arc_curve_radi = 0.01,
      face_normal_vertex_gm_radi = 0.03,
      gm_distance = 2.1,
      gm_radi = 1.;
int arcs_seg_count = 13;

double ground_box_y = -3;
Vector3 ground_box_shape({10,1,10});

polyscope::SlicePlane* psPlane;


std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr;
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;


polyscope::SurfaceMesh *psMesh, *dummy_psMesh2;
polyscope::PointCloud *gauss_map_pc, *face_normals_pc;


void update_positions(){
    VertexData<Vector3> new_positions = my_env->get_new_positions();
    psMesh->updateVertexPositions(new_positions.raw());
}


void initialize_vis(){
    psMesh = polyscope::registerSurfaceMesh("my polyhedra", geometry->inputVertexPositions, mesh->getFaceVertexList());
    // ground plane on Polyscope has a weird height setting (scaled factor..)
    psPlane = polyscope::addSceneSlicePlane();
    psPlane->setDrawPlane(true);  // render the semi-transparent gridded plane
    psPlane->setDrawWidget(false);
    psPlane->setPose(glm::vec3{0., ground_box_y + 1, 0.}, glm::vec3{0., 1., 0.});
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
  // just draw the sphere next to the main surface
  Vector3 shift = {0., 0. , gm_distance};
  std::vector<Vector3> sphere_pos = {shift};
  gauss_map_pc = polyscope::registerPointCloud("Gauss Map", sphere_pos);
  gauss_map_pc->setPointColor({0.74,0.7,0.9});
  gauss_map_pc->setPointRadius(gm_radi, false);
  gauss_map_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

  // point cloud for face normals
  std::vector<Vector3> face_normal_points;
  for (Face f: mesh->faces()){
    Vector3 normal_pos_on_gm = -geometry->faceNormal(f) + shift;
    face_normal_points.push_back(normal_pos_on_gm);
  }
  face_normals_pc = polyscope::registerPointCloud("Face Normals", face_normal_points);
  face_normals_pc->setPointRadius(face_normal_vertex_gm_radi, false);
  face_normals_pc->setPointColor({0.9,0.9,0.9});
  face_normals_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  
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
    Vector3 n1 = -geometry->faceNormal(f1),
            n2 = -geometry->faceNormal(f2);
    // draw with polyscope
    draw_arc_on_sphere(n1, n2, shift, gm_radi, arcs_seg_count, e.getIndex());
  }
}




// polyscope callback
void myCallback() {
    // if (ImGui::Button("Take simulation steps")){
    //     take_simulation_steps();
    // }
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

    if (ImGui::Button("refresh")){
        my_env->refresh(find_center_of_mass(*mesh, *geometry), Vector3({0,0,1}));
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
}


int main(int argc, char* argv[])
{
#ifdef BT_USE_DOUBLE_PRECISION
	printf("BT_USE_DOUBLE_PRECISION\n");
#else
    printf("Single precision\n");
#endif
    std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra("tet2");
    mesh = mesh_ptr.release(); 
    geometry = geometry_ptr.release();
    // center_and_normalize(mesh, geometry);

    my_env = new PhysicsEnv();
	my_env->init_physics();
    my_env->init_geometry(mesh, geometry);

	///-----initialization_end-----
    my_env->add_ground(ground_box_y, ground_box_shape);
    my_env->add_object(find_center_of_mass(*mesh, *geometry), Vector3({1,0,0}));

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
