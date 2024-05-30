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
#pragma once

#include "LinearMath/btVector3.h"
#include "btBulletDynamicsCommon.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;


geometrycentral::DenseMatrix<double> openGL_mat_to_GC_mat(btScalar* m);
btQuaternion quaternion_from_g_vec(Vector3 default_g_vec, Vector3 g_vec);

Vector<Vector3> apply_trans_to_positions(Vector<Vector3> init_positions, geometrycentral::DenseMatrix<double> trans_mat);

class PhysicsEnv {
  public:
    // GC stuff
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;
    
    // keeps track of current positions via a transformation
    Vector3 initial_orientation;
    btTransform init_btTransform;
    btTransform current_btTrans;
    std::vector<Vector3> orientation_trail;
    std::vector<geometrycentral::DenseMatrix<double>> trans_mat_trail;
    
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

    int MAX_ITERS = 100000, MIN_ITERS = 100;
    double tol = 1e-9,
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
    void add_object(Vector3 center_of_mass, Vector3 orientation);
    // just delete the polyhedra
    void delete_object();

    // take a time step in the simulation
    void take_step(int step_count, double step_size);

    // find the bottom face
    Face get_touching_face(VertexData<Vector3> new_positions);
    // simulate until reaching the stable orientation
    Face final_stable_face(bool save_trail = false);


    //restart with the same polyhedra
    void refresh(Vector3 G, Vector3 g_vec);
    Vector3 get_current_G();
    Vector3 get_current_orientation();
    // delete stuff; for cache reasons?
    void delete_bullet_objects();

    // 
    // get btTransform object as GC matrix
    geometrycentral::DenseMatrix<double> btTrans_to_GC_Mat(btTransform current_btTrans);
    Vector<Vector3> get_new_positions(Vector<Vector3> old_poses);
};