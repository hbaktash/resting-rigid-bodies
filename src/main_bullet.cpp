#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "LinearMath/btVector3.h"
#include "btBulletDynamicsCommon.h"
// #include <stdio.h>

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

size_t MAX_ITERS = 200, PER_STEP_ITERS = 1, iter = 0;
Vector3 init_position({0.,0.,0.});
Vector3 ground_box_shape({1,1,1}), 
        ground_origin({0,-4,0});

polyscope::SlicePlane* psPlane;


std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr;
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

polyscope::SurfaceMesh *psMesh;

// Bulet3 vars
btConvexHullShape* bt_my_polyhedra;

btDefaultCollisionConfiguration* collisionConfiguration;
btCollisionDispatcher* dispatcher;
btBroadphaseInterface* overlappingPairCache;
btSequentialImpulseConstraintSolver* solver;
btDiscreteDynamicsWorld* dynamicsWorld;
//keep track of the shapes, we release memory at exit.
//make sure to re-use collision shapes among rigid bodies whenever possible!
btAlignedObjectArray<btCollisionShape*> collisionShapes;

// for future display
std::vector<Vector3> sphere_positions;


geometrycentral::DenseMatrix<double> openGL_mat_to_GC_mat(btScalar* m){
    geometrycentral::DenseMatrix<double> trans_mat(4,4);
    for (int flat_ind = 0; flat_ind < 16; flat_ind++){
        double tmp_var = float(m[flat_ind]);
        int i = flat_ind/4, j = flat_ind % 4;
        trans_mat.coeffRef(j, i) = tmp_var;
    }
    std::cout<< trans_mat <<"\n";
    return trans_mat.inverse();
}


void take_simulation_step(){
    //print positions of all objects
    for (int iter = 0; iter < PER_STEP_ITERS; iter++){
        dynamicsWorld->stepSimulation(1./30., 10);
        for (int j = dynamicsWorld->getNumCollisionObjects() - 1; j >= 0; j--) {
            btCollisionObject* obj = dynamicsWorld->getCollisionObjectArray()[j];
            btRigidBody* body = btRigidBody::upcast(obj);
            btTransform trans;
            if (body && body->getMotionState())
            {
                body->getMotionState()->getWorldTransform(trans);
            }
            else
            {
                trans = obj->getWorldTransform();
            }
            if (j == 1){// sphere
                btScalar m[16];
                trans.getOpenGLMatrix(m);
                geometrycentral::DenseMatrix<double> trans_mat = openGL_mat_to_GC_mat(m);
                // btVector3 rot = trans.getRotation();
                // btConvexHullShape* tmp_convex_hull = (btConvexHullShape*) obj;
                // const btVector3* points = bt_my_polyhedra->getPoints();
                // std::cout<<points->getX()<<"\n";
                std::vector<Vector3> new_positions;
                for (Vertex v: mesh->vertices()) {
                    Vector3 old_p = geometry->inputVertexPositions[v];
                    Vector<double> tmp_vec(4);
                    tmp_vec[0] = old_p.x; tmp_vec[1] = old_p.y; tmp_vec[2] = old_p.z; tmp_vec[3] = 1.;
                    tmp_vec = trans_mat*tmp_vec;
                    Vector3 new_p({tmp_vec[0]/tmp_vec[3], tmp_vec[1]/tmp_vec[3], tmp_vec[2]/tmp_vec[3]});
                    new_positions.push_back(new_p);
                }
                psMesh->updateVertexPositions(new_positions);
                // printf("world pos object %d = %f,%f,%f\n", j, float(trans.getOrigin().getX()), float(trans.getOrigin().getY()), float(trans.getOrigin().getZ()));
                // Vector3 curr_vec = Vector3({float(trans.getOrigin().getX()), float(trans.getOrigin().getY()), float(trans.getOrigin().getZ())});
                // std::vector<Vector3> curr_posz = {curr_vec};
                // sphere_positions.push_back(curr_vec);
                // psSphere = polyscope::registerPointCloud("Dummy sphere", curr_posz);
                // psSphere->updatePointPositions(curr_posz);
            }
        }
    }
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
    // ground plane on Polyscope
    // psPlane = polyscope::addSceneSlicePlane();
    // psPlane->setDrawPlane(true);  // render the semi-transparent gridded plane
    // psPlane->setDrawWidget(true);
    // psPlane->setPose(glm::vec3{0., -3., 0.}, glm::vec3{0., -1., 0.});
    // aAdd to sphere polyscope 
    // std::vector<Vector3> sphere_pos = {init_position};
    // psSphere = polyscope::registerPointCloud("Dummy sphere", sphere_pos);
    // psSphere->setPointColor({0.4, 0.2, 0.2});
    // psSphere->setPointRadius(sphere_radi, true);
    // psSphere->setPointRenderMode(polyscope::PointRenderMode::Sphere);
    // psSphere->setCullWholeElements(false);
    psMesh = polyscope::registerSurfaceMesh("my polyhedra", geometry->inputVertexPositions, mesh->getFaceVertexList());
}


// polyscope callback
void myCallback() {
    // if (ImGui::Button("Take simulation steps")){
    //     take_simulation_steps();
    // }
    if (ImGui::Button("take simulation step")){
        take_simulation_step();
    }
    if (ImGui::Button("Clean up stuff")){
        delete_stuff();
    }
}


void initialize_bullet_vars(){
    ///-----initialization_start-----
	///collision configuration contains default setup for memory, collision setup. Advanced users can create their own configuration.
	collisionConfiguration = new btDefaultCollisionConfiguration();
	///use the default collision dispatcher. For parallel processing you can use a diffent dispatcher (see Extras/BulletMultiThreaded)
	dispatcher = new btCollisionDispatcher(collisionConfiguration);
	///btDbvtBroadphase is a good general purpose broadphase. You can also try out btAxis3Sweep.
	overlappingPairCache = new btDbvtBroadphase();
	///the default constraint solver. For parallel processing you can use a different solver (see Extras/BulletMultiThreaded)
	solver = new btSequentialImpulseConstraintSolver;
	dynamicsWorld = new btDiscreteDynamicsWorld(dispatcher, overlappingPairCache, solver, collisionConfiguration);
}

//the ground is a cube of side 100 at position y = -1.
//the sphere will hit it at y = -6, with center at -5	
void add_ground(){
    btCollisionShape* groundShape = new btBoxShape(btVector3(btScalar(1.), btScalar(1.), btScalar(1.)));

    collisionShapes.push_back(groundShape);

    btTransform groundTransform;
    groundTransform.setIdentity();
    groundTransform.setOrigin(btVector3(0, -1., 0));

    btScalar mass(0.);

    //rigidbody is dynamic if and only if mass is non zero, otherwise static
    bool isDynamic = (mass != 0.f);

    btVector3 localInertia(0, 0, 0);
    if (isDynamic)
        groundShape->calculateLocalInertia(mass, localInertia);

    //using motionstate is optional, it provides interpolation capabilities, and only synchronizes 'active' objects
    btDefaultMotionState* myMotionState = new btDefaultMotionState(groundTransform);
    btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, groundShape, localInertia);
    btRigidBody* body = new btRigidBody(rbInfo);

    //add the body to the dynamics world
    dynamicsWorld->addRigidBody(body);
    
}


void add_object(){
    //create a dynamic rigidbody
    //btCollisionShape* colShape = new btBoxShape(btVector3(1,1,1));
    std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra("tet");
    mesh = mesh_ptr.release(); geometry = geometry_ptr.release();
    
    bt_my_polyhedra = new btConvexHullShape();
    for (Vertex v: mesh->vertices()) {
        Vector3 p = geometry->inputVertexPositions[v];
        btVector3 btP(p.x, p.y, p.z);
        bt_my_polyhedra->addPoint(btP);
    }
    btCollisionShape* colShape = new btConvexHullShape(*bt_my_polyhedra);
    
    // btCollisionShape* colShape = new btSphereShape(btScalar(sphere_radi));
    collisionShapes.push_back(colShape);

    /// Create Dynamic Objects
    btTransform startTransform;
    startTransform.setIdentity();

    btScalar mass(1.f);

    //rigidbody is dynamic if and only if mass is non zero, otherwise static
    bool isDynamic = (mass != 0.f);

    btVector3 localInertia(0, 0, 0);
    if (isDynamic)
        colShape->calculateLocalInertia(mass, localInertia);

    startTransform.setOrigin(btVector3(init_position.x, init_position.y, init_position.z));

    //using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
    btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
    btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, colShape, localInertia);
    btRigidBody* body = new btRigidBody(rbInfo);

    dynamicsWorld->addRigidBody(body);
}


int main(int argc, char* argv[])
{
#ifdef BT_USE_DOUBLE_PRECISION
	printf("BT_USE_DOUBLE_PRECISION\n");
#else
        printf("Single precision\n");
#endif

	initialize_bullet_vars();
	dynamicsWorld->setGravity(btVector3(0, -10, 0));

	///-----initialization_end-----
    add_ground();
    add_object();

	///-----stepsimulation_end-----
    // take_simulation_step();

    polyscope::init();
    // polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    // polyscope::view::upDir = polyscope::view::UpDir::NegYUp;
    // polyscope::options::groundPlaneHeightFactor = -3.; // adjust the plane height
    polyscope::state::userCallback = myCallback;
    initialize_vis();
    
    // Set the callback function


    // Give control to the polyscope gui
    polyscope::show();

    
    return EXIT_SUCCESS;
	
}
