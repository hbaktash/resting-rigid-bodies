#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
// #include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
// #include "bullet3/examples/BasicExample.h"
#include "args/args.hxx"
#include "imgui.h"

#include "mesh_factory.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

double sphere_radi = 0.2;
Vector3 init_position({0., 4., 0.}), other_pos({0., 3.1, 0.});
bool other = false;
polyscope::PointCloud *psSphere;
polyscope::SlicePlane* psPlane;

std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr;
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

polyscope::SurfaceMesh *psMesh;

void initialize_vis(){
    // ground plane on Polyscope
    // psPlane = polyscope::addSceneSlicePlane();
    // psPlane->setDrawPlane(true);  // render the semi-transparent gridded plane
    // psPlane->setDrawWidget(true);
    // psPlane->setPose(glm::vec3{0., -1.1, 0.}, glm::vec3{0., -1., 0.});
    // aAdd to sphere polyscope 
    std::vector<Vector3> sphere_pos = {init_position};
    psSphere = polyscope::registerPointCloud("Dummy sphere", sphere_pos);
    psSphere->setPointColor({0.4, 0.2, 0.2});
    psSphere->setPointRadius(sphere_radi, true);
    psSphere->setPointRenderMode(polyscope::PointRenderMode::Sphere);
    psSphere->setCullWholeElements(false);

    psMesh = polyscope::registerSurfaceMesh("some mesh", geometry->inputVertexPositions, mesh->getFaceVertexList());
}


void move_sphere(){
    Vector3 curr_vec = other ? init_position : other_pos;
    other = !other;
    std::vector<Vector3> curr_posz = {curr_vec};
    psSphere->updatePointPositions(curr_posz);


}

// polyscope callback
void myCallback() {
    // if (ImGui::Button("Take simulation steps")){
    //     take_simulation_steps();
    // }
    if (ImGui::Button("Move sphere")){
        move_sphere();
    }
}


int main(int argc, char* argv[]) {

    polyscope::init();
    // polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    // polyscope::view::upDir = polyscope::view::UpDir::NegYUp;
    // polyscope::options::groundPlaneHeightFactor = -3.; // adjust the plane height
    polyscope::state::userCallback = myCallback;
    std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra("tet");
    mesh = mesh_ptr.release(); geometry = geometry_ptr.release();
    initialize_vis();
    
    // Set the callback function


    // Give control to the polyscope gui
    polyscope::show();

    
    return EXIT_SUCCESS;
	
}
