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


PhysicsEnv* my_env;
float step_size = 0.03;
int step_count = 1;

double ground_box_y = -3;
Vector3 ground_box_shape({10,1,10});

polyscope::SlicePlane* psPlane;


std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr;
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

polyscope::SurfaceMesh *psMesh;


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
    
    if (ImGui::Button("update positions")){
        update_positions();
    }

    if (ImGui::Button("refresh")){
        my_env->refresh(find_center_of_mass(*mesh, *geometry), Vector3({0,0,1}));
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
    center_and_normalize(mesh, geometry);

    my_env = new PhysicsEnv();
	my_env->init_physics();
    my_env->init_geometry(mesh, geometry);

	///-----initialization_end-----
    my_env->add_ground(ground_box_y, ground_box_shape);
    my_env->add_object(find_center_of_mass(*mesh, *geometry), Vector3({1,0,0}));

    Face touching_face = my_env->final_stable_face(Vector3({0,-1,0}));
    printf("final touching face is %d\n", touching_face.getIndex());
    std::cout << "face normal is "<< geometry->faceNormal(touching_face)<< "\n";

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
