#pragma once

#include "LinearMath/btVector3.h"
#include "btBulletDynamicsCommon.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;


geometrycentral::DenseMatrix<double> openGL_mat_to_GC_mat(btScalar* m);
btQuaternion quaternion_from_g_vec(Vector3 default_g_vec, Vector3 g_vec);

class PhysicsEnv {
  public:
    // GC stuff
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
    
    // keeps track of current positions via a transformation
    btTransform current_btTrans;

    // Bullet3 shape var
    btConvexHullShape* bt_my_polyhedra;

    // Bullet3 env; collisions, dynamic objects, ..
    btDefaultCollisionConfiguration* collisionConfiguration;
    btCollisionDispatcher* dispatcher;
    btBroadphaseInterface* overlappingPairCache;
    btSequentialImpulseConstraintSolver* solver;
    btDiscreteDynamicsWorld* dynamicsWorld;
    //keep track of the shapes, we release memory at exit.
    //make sure to re-use collision shapes among rigid bodies whenever possible!
    btAlignedObjectArray<btCollisionShape*> collisionShapes;

    int MAX_ITERS = 200;
    double tol = 1e-8,
           default_step_size = 0.017; // 1/60


    // Constructor
    PhysicsEnv(){}
    // Destructor
    ~PhysicsEnv(){ 
      delete_bullet_objects(); 
    }

    // initialize physics and bullet vars, and mesh+geometry
    void init_physics();
    void init_geometry(ManifoldSurfaceMesh* convex_mesh, VertexPositionGeometry* convex_geometry);

    // add static ground
    void add_ground(double ground_box_y, Vector3 ground_box_shape);

    // add convex object; TODO double check center of mass issue
    void add_object(Vector3 center_of_mass, Vector3 g_vec);

    // take a time step in the simulation
    void take_step(int step_count, double step_size);

    // simulate until reaching the stable orientation
    Face final_stable_face(Vector3 g_vec);

    // delete stuff; for cache reasons?
    void delete_bullet_objects();


    // get btTransform object as GC matrix
    std::vector<Vector3> get_new_positions();
};