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

void PhysicsEnv::take_step(double step_size){
    //print positions of all objects
    dynamicsWorld->stepSimulation(step_size, 10);
    for (int j = dynamicsWorld->getNumCollisionObjects() - 1; j >= 0; j--) {
        btCollisionObject* obj = dynamicsWorld->getCollisionObjectArray()[j];
        btRigidBody* body = btRigidBody::upcast(obj);
        btTransform trans;
        if (body && body->getMotionState()) {
            body->getMotionState()->getWorldTransform(trans);
        }
        else {
            trans = obj->getWorldTransform();
        }
        if (j == 1){ // object of interest; i.e. not ground
            btVector3 center_of_mass = body->getCenterOfMassPosition();
            // std::cout<< "center of mass "<<center_of_mass.getX()<<","<< center_of_mass.getY()<< "," << center_of_mass.getZ() << "\n";
            current_btTrans = trans;
            
            
            // printf("world pos object %d = %f,%f,%f\n", j, float(trans.getOrigin().getX()), float(trans.getOrigin().getY()), float(trans.getOrigin().getZ()));
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


void PhysicsEnv::add_ground(){
    btCollisionShape* groundShape = new btBoxShape(btVector3(btScalar(10.), btScalar(1.), btScalar(10.)));
    collisionShapes.push_back(groundShape);
    btTransform groundTransform;
    groundTransform.setIdentity();
    groundTransform.setOrigin(btVector3(0, -3., 0));
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


void PhysicsEnv::add_object(ManifoldSurfaceMesh* convex_mesh, VertexPositionGeometry* convex_geometr, Vector3 G){
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
