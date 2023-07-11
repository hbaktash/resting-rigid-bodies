#include "bullet_sim.h"
#include "geometry_utils.h"

geometrycentral::DenseMatrix<double> openGL_mat_to_GC_mat(btScalar* m){
    geometrycentral::DenseMatrix<double> trans_mat(4,4);
    for (int flat_ind = 0; flat_ind < 16; flat_ind++){
        double tmp_var = float(m[flat_ind]);
        int i = flat_ind/4, j = flat_ind % 4;
        trans_mat.coeffRef(i, j) = tmp_var;
    }
    // std::cout<< trans_mat <<"\n";
    return trans_mat.transpose();
}

btQuaternion quaternion_from_g_vec(Vector3 default_g_vec, Vector3 g_vec){
    Vector3 plane_normal = cross(g_vec, default_g_vec).normalize();
    double angle = acos(dot(g_vec, default_g_vec)/(norm(g_vec)*norm(default_g_vec)));
    return btQuaternion(btVector3(plane_normal.x, plane_normal.y, plane_normal.z), btScalar(angle));
}


std::vector<Vector3> PhysicsEnv::get_new_positions(){
    btScalar m[16];
    current_btTrans.getOpenGLMatrix(m);
    geometrycentral::DenseMatrix<double> trans_mat = openGL_mat_to_GC_mat(m);
    std::vector<Vector3> new_positions;
    for (Vertex v: mesh->vertices()) {
        Vector3 old_p = geometry->inputVertexPositions[v];
        Vector<double> tmp_vec(4);
        tmp_vec[0] = old_p.x; tmp_vec[1] = old_p.y; tmp_vec[2] = old_p.z; tmp_vec[3] = 1.;
        tmp_vec = trans_mat*tmp_vec;
        Vector3 new_p({tmp_vec[0]/tmp_vec[3], tmp_vec[1]/tmp_vec[3], tmp_vec[2]/tmp_vec[3]});
        new_positions.push_back(new_p);
    }
    return new_positions;
}

void PhysicsEnv::take_step(int step_count, double step_size){
    //print positions of all objects
    for (int i = 0; i < step_count; i++){
        dynamicsWorld->stepSimulation(step_size, 10);
    }
    for (int j = dynamicsWorld->getNumCollisionObjects() - 1; j >= 0; j--) {
        if (j == 1){ // object of interest; i.e. not ground
            btCollisionObject* obj = dynamicsWorld->getCollisionObjectArray()[j];
            btRigidBody* body = btRigidBody::upcast(obj);
            btTransform trans;
            if (body && body->getMotionState()) {
                body->getMotionState()->getWorldTransform(trans);
            }
            else {
                trans = obj->getWorldTransform();
            }
            // btVector3 center_of_mass = body->getCenterOfMassPosition();
            // std::cout<< "center of mass "<<center_of_mass.getX()<<","<< center_of_mass.getY()<< "," << center_of_mass.getZ() << "\n";
            current_btTrans = trans;
        }
    }
}

void PhysicsEnv::init_physics(){
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

    // probably do this at construction time? maybe this whole function?
    dynamicsWorld->setGravity(btVector3(0, -10, 0));
}


void PhysicsEnv::init_geometry(ManifoldSurfaceMesh* convex_mesh, VertexPositionGeometry* convex_geometry){
    this->mesh = convex_mesh;
    this->geometry = convex_geometry;
}

void PhysicsEnv::add_ground(double ground_box_y, Vector3 ground_box_shape){
    btCollisionShape* groundShape = new btBoxShape(btVector3(ground_box_shape.x, ground_box_shape.y, ground_box_shape.z));
    collisionShapes.push_back(groundShape);
    btTransform groundTransform;
    groundTransform.setIdentity();
    groundTransform.setOrigin(btVector3(0, ground_box_y, 0));
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


void PhysicsEnv::add_object(Vector3 G, Vector3 g_vec){
    //create a dynamic rigidbody
    //btCollisionShape* colShape = new btBoxShape(btVector3(1,1,1));
    bt_my_polyhedra = new btConvexHullShape();
    for (Vertex v: mesh->vertices()) {
        Vector3 p = geometry->inputVertexPositions[v];
        btVector3 btP(p.x, p.y, p.z);
        bt_my_polyhedra->addPoint(btP);
    }
    btCollisionShape* colShape = new btConvexHullShape(*bt_my_polyhedra);
    collisionShapes.push_back(colShape);

    /// Create Dynamic Objects
    btTransform startTransform;
    startTransform.setIdentity();
    // rotate to face the g_vec down
    startTransform.setRotation(quaternion_from_g_vec(Vector3({0., -1., 0.}), g_vec));
    // any non-zero mass works
    btScalar mass(1.f);

    //rigidbody is dynamic if and only if mass is non zero, otherwise static
    bool isDynamic = (mass != 0.f);

    btVector3 localInertia(0, 0, 0);
    if (isDynamic)
        colShape->calculateLocalInertia(mass, localInertia);
    startTransform.setOrigin(btVector3(G.x, G.y, G.z));

    //using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
    btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
    btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, colShape, localInertia);
    btRigidBody* body = new btRigidBody(rbInfo);

    dynamicsWorld->addRigidBody(body);
}


// take many steps till stable
Face PhysicsEnv::final_stable_face(Vector3 initial_g_vec){
    Vector3 old_G, current_G;
    
}


void PhysicsEnv::delete_bullet_objects(){
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
    //delete containers
    delete dynamicsWorld;
    delete solver;
    delete overlappingPairCache;
    delete dispatcher;
    delete collisionConfiguration;
    //next line is optional: it will be cleared by the destructor when the array goes out of scope
    collisionShapes.clear();
}
