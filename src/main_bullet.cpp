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

size_t MAX_ITERS = 200, PER_STEP_ITERS = 1, iter = 0;
Vector3 init_position({0.,0.,0.}),
        G({0.,0.,0.});
Vector3 ground_box_shape({1,1,1}), 
        ground_origin({0,-7,0});

polyscope::SlicePlane* psPlane;


std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr;
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

polyscope::SurfaceMesh *psMesh;

// // Bulet3 vars
// btConvexHullShape* bt_my_polyhedra;

// btDefaultCollisionConfiguration* collisionConfiguration;
// btCollisionDispatcher* dispatcher;
// btBroadphaseInterface* overlappingPairCache;
// btSequentialImpulseConstraintSolver* solver;
// btDiscreteDynamicsWorld* dynamicsWorld;
// //keep track of the shapes, we release memory at exit.
// //make sure to re-use collision shapes among rigid bodies whenever possible!
// btAlignedObjectArray<btCollisionShape*> collisionShapes;


// geometrycentral::DenseMatrix<double> openGL_mat_to_GC_mat(btScalar* m){
//     geometrycentral::DenseMatrix<double> trans_mat(4,4);
//     for (int flat_ind = 0; flat_ind < 16; flat_ind++){
//         double tmp_var = float(m[flat_ind]);
//         int i = flat_ind/4, j = flat_ind % 4;
//         trans_mat.coeffRef(i, j) = tmp_var;
//     }
//     // std::cout<< trans_mat <<"\n";
//     return trans_mat.transpose();
// }


void update_positions(){
    std::vector<Vector3> new_positions = my_env->get_new_positions();
    psMesh->updateVertexPositions(new_positions);
}


void delete_stuff(){
    //cleanup in the reverse order of creation/initialization
	///-----cleanup_start-----
	//remove the rigidbodies from the dynamics world and delete them
    for (int i = dynamicsWorld->getNumCollisionObjects() - 1; i >= 0; i--) {
        btCollisionObject* obj = dynamicsWorld->getCollisionObjectArray()[i];
        btRigidBody* body = btRigidBody::upcast(obj);
        if (body && body->getMotionState())
        {
            delete body->getMotionState();
        }
        dynamicsWorld->removeCollisionObject(obj);
        delete obj;
    }
    //delete collision shapes
    for (int j = 0; j < collisionShapes.size(); j++) {
        btCollisionShape* shape = collisionShapes[j];
        collisionShapes[j] = 0;
        delete shape;
    }

    //delete dynamics world
    delete dynamicsWorld;

    //delete solver
    delete solver;

    //delete broadphase
    delete overlappingPairCache;

    //delete dispatcher
    delete dispatcher;

    delete collisionConfiguration;

    //next line is optional: it will be cleared by the destructor when the array goes out of scope
    collisionShapes.clear();
}

void initialize_vis(){

    psMesh = polyscope::registerSurfaceMesh("my polyhedra", geometry->inputVertexPositions, mesh->getFaceVertexList());
    
    // ground plane on Polyscope?
    // psPlane = polyscope::addSceneSlicePlane();
    // psPlane->setDrawPlane(true);  // render the semi-transparent gridded plane
    // psPlane->setDrawWidget(true);
    // psPlane->setPose(glm::vec3{0., -3., 0.}, glm::vec3{0., -1., 0.});
}


// polyscope callback
void myCallback() {
    // if (ImGui::Button("Take simulation steps")){
    //     take_simulation_steps();
    // }
    if (ImGui::Button("update positions")){
        update_positions();
    }
    if (ImGui::Button("Clean up stuff")){
        delete_stuff();
    }
}


// ground is a cube
void add_ground(){
    
    
}


void add_object(){
    std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra("tet2");
    mesh = mesh_ptr.release(); geometry = geometry_ptr.release();
}


int main(int argc, char* argv[])
{
#ifdef BT_USE_DOUBLE_PRECISION
	printf("BT_USE_DOUBLE_PRECISION\n");
#else
        printf("Single precision\n");
#endif
    my_env = new PhysicsEnv();
	my_env->init_physics();
	

	///-----initialization_end-----
    my_env->add_ground();
    my_env->add_object(mesh, geometry, find_center_of_mass(*mesh, *geometry));

	///-----stepsimulation_end-----
    // take_simulation_step();

    polyscope::init();
    // polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    // polyscope::view::upDir = polyscope::view::UpDir::YUp;
    // polyscope::options::groundPlaneHeightFactor = -3.; // adjust the plane height
    polyscope::state::userCallback = myCallback;
    initialize_vis();
    
    // Set the callback function


    // Give control to the polyscope gui
    polyscope::show();

    
    return EXIT_SUCCESS;
	
}
