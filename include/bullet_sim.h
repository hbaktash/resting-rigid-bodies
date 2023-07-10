#pragma once

#include "LinearMath/btVector3.h"
#include "btBulletDynamicsCommon.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;


geometrycentral::DenseMatrix<double> openGL_mat_to_GC_mat(btScalar* m);


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

    // initialize physics and bullet vars
    void init_physics();

    // add static ground
    void add_ground();

    // add convex object; TODO double check center of mass issue
    void add_object(ManifoldSurfaceMesh* convex_mesh, VertexPositionGeometry* convex_geometry, Vector3 center_of_mass);

    // take a time step in the simulation
    void take_step(double step_size);

    // delete stuff; for cache reasons?
    void delete_bullet_objects();


    // get btTransform object as GC matrix
    std::vector<Vector3> get_new_positions();
};