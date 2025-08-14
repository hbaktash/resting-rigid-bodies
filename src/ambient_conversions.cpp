#include "ambient_conversions.h"


std::pair<std::vector<Eigen::Matrix4d>, std::vector<Vector3>> 
generate_transformations_for_orientation_sequence(Vector3 initial_orientation,
    Forward3DSolver* forwardSolver, Vector3 floor_vec, double goal_angle_step){
    // Get the snail trail and make forward log
    std::vector<Vector3> snail_trail = forwardSolver->quasi_static_drop(initial_orientation);
    // split the snail trail
    std::vector<Vector3> snail_trail_refined;
    for (int i = 1; i < snail_trail.size()-1; i++){
        Vector3 local_axis = cross(snail_trail[i-1], snail_trail[i]).normalize();
        double local_total_angle = angle(snail_trail[i-1], snail_trail[i]);
        int steps = (int)ceil(local_total_angle/goal_angle_step); // + 1
        // in steps = 2;
        for (int t = 0; t < steps; t++){
            double angle_0 = local_total_angle * (double)t/double(steps);
            Vector3 normal_0 = snail_trail[i-1].rotateAround(local_axis, angle_0);
            snail_trail_refined.push_back(normal_0);
        }
    }
    snail_trail_refined.push_back(snail_trail[snail_trail.size()-1]);

    VertexData<Vector3> init_hull_positions = forwardSolver->hullGeometry->inputVertexPositions;
    VertexData<Vector3> tmp_hull_positions(*forwardSolver->hullMesh);
    tmp_hull_positions = init_hull_positions;
    
    std::vector<Vector3> saved_snail_trail_refined = snail_trail_refined; // for future output
    // get transformations
    std::vector<Eigen::Matrix4d> transformations;
    // transformations.push_back(Eigen::Matrix4d::Identity()); // initial identity matrix
    for (int i = 0; i < snail_trail_refined.size(); i++){
        // get n_i locally
        // SO3 conversion
        Vector3 normal_0 = snail_trail_refined[i];
        Vector3 rot_axis = cross(normal_0, floor_vec).normalize();
        double rot_angle = angle(normal_0, floor_vec);

        for (int j = 0; j < snail_trail_refined.size(); j++){
            snail_trail_refined[j] = snail_trail_refined[j].rotateAround(rot_axis, rot_angle);
        }
        // shift contact to origin
        double lowest_height = 1e8;
        Vector3 contact_p;
        for (Vertex v: forwardSolver->hullMesh->vertices()){
            double tmp_height = dot(tmp_hull_positions[v], -floor_vec);
            if (tmp_height < lowest_height){
                lowest_height = tmp_height;
                contact_p = tmp_hull_positions[v];
            }
        }
        Eigen::AngleAxisd aa(rot_angle, Eigen::Vector3d(rot_axis.x, rot_axis.y, rot_axis.z));
        Eigen::Matrix3d rotation_matrix = aa.toRotationMatrix();
        // do the rotation around the contact point
        // shifting contact
        tmp_hull_positions -= contact_p;
        for (Vertex v: forwardSolver->hullMesh->vertices()){
            // tmp_hull_positions[v] = tmp_hull_positions[v].rotateAround(rot_axis, rot_angle);
            Eigen::Vector3d tmp_v(tmp_hull_positions[v].x, tmp_hull_positions[v].y, tmp_hull_positions[v].z);
            tmp_v = rotation_matrix * tmp_v;
            tmp_hull_positions[v] = Vector3({tmp_v[0], tmp_v[1], tmp_v[2]});
        }
        // shift back
        tmp_hull_positions += contact_p;
        // axis angle rotation to matrix

        // TODO: this is a hack, should be fixed in the future
        // height shift; not necessary if the contact point is chosen correctly; tricky when on an unstable triangle
        lowest_height = 1e8;
        for (Vertex v: forwardSolver->hullMesh->vertices()){
            if (dot(tmp_hull_positions[v], -floor_vec) < lowest_height) //tmp_hull_positions[v].y < lowest_height
                lowest_height = dot(tmp_hull_positions[v], -floor_vec);
        }
        Eigen::Matrix4d transformation_matrix;
        transformation_matrix.setIdentity();
        transformation_matrix.block<3,3>(0,0) = rotation_matrix;
        // the last column is the translation
        Eigen::Vector3d contact_p_eigen(contact_p.x, contact_p.y, contact_p.z),
                      floor_vec_eigen(floor_vec.x, floor_vec.y, floor_vec.z);
        Eigen::Vector3d translation = -rotation_matrix * contact_p_eigen + 
                                      contact_p_eigen + lowest_height * floor_vec_eigen; // + Eigen::Vector3d(0, -lowest_height, 0)
        transformation_matrix.block<3,1>(0,3) = translation;
        transformations.push_back(transformation_matrix);
        // update tmp_hull_positions height
        for (Vertex v: forwardSolver->hullMesh->vertices()){
            tmp_hull_positions[v] += floor_vec * lowest_height; //tmp_hull_positions[v].y += -lowest_height
        }
    }
    return {transformations, saved_snail_trail_refined};
}


std::pair<std::vector<Eigen::Matrix4d> , std::vector<Vector3>>
QS_trail_to_global_transformations(
    Forward3DSolver* forwardSolver, Vector3 initial_orientation, 
    Vector3 down_vec, double min_angle_step
){
    // initialize the trail in forwardSolver
    std::vector<Vector3> snail_trail = forwardSolver->quasi_static_drop(initial_orientation);
    // split the snail trail
    auto [local_transformation_matrices, saved_snail_trail_refined] = 
        generate_transformations_for_orientation_sequence(
        initial_orientation, forwardSolver, down_vec, min_angle_step);
    // local to global; L_n = R_n * R_(n-1) * ... * R_0
    std::vector<Eigen::Matrix4d> transformation_matrices;
    Eigen::Matrix4d global_transformation = Eigen::Matrix4d::Identity();
    for (int i = 0; i < local_transformation_matrices.size(); i++){
        Eigen::Matrix4d local_transformation = local_transformation_matrices[i];
        global_transformation = local_transformation * global_transformation;
        transformation_matrices.push_back(global_transformation);
    }
    
    return {transformation_matrices, saved_snail_trail_refined};
}


VertexData<Vector3> orient_mesh_to_unit_vec(Vector3 orientation, Vector3 down_vec,
                               ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry){
    Vector3 axis = cross(orientation, down_vec);
    if (axis.norm() < 1e-6) // parallel
        axis = Vector3({1,0,0});
    axis = axis.normalize();
    double angle = acos(dot(orientation, down_vec));
    VertexData<Vector3> new_positions(*mesh);
    for (Vertex v: mesh->vertices()){
        Vector3 p = geometry->inputVertexPositions[v];
        new_positions[v] = p.rotateAround(axis, angle);
    }
    return new_positions;
}


// get trans mat for orientating mesh to unit vector
Eigen::Matrix4d trans_mat_for_orientation(Vector3 orientation, Vector3 down_vec,
                                          ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry){
    Vector3 axis = cross(orientation, down_vec);
    if (axis.norm() < 1e-6) // parallel
        axis = Vector3({1,0,0});
    axis = axis.normalize();
    double angle = acos(dot(orientation, down_vec));
    Eigen::AngleAxisd aa(angle, Eigen::Vector3d(axis.x, axis.y, axis.z));
    Eigen::Matrix3d rotation_matrix = aa.toRotationMatrix();
    Eigen::Matrix4d transformation_matrix;
    transformation_matrix.setIdentity();
    transformation_matrix.block<3,3>(0,0) = rotation_matrix;
    // find translation for the contact point to touch the ground
    double max_dist = -1e8;
    Vector3 contact_point;
    for (Vertex v: mesh->vertices()){
        Vector3 p = geometry->inputVertexPositions[v];
        double dist = dot(p, orientation);
        if (dist > max_dist){
            max_dist = dist;
            contact_point = p;
        }
    }
    // find the contact point
    Vector3 shift_vec = max_dist * down_vec;
    Eigen::Vector3d shift_eigen(shift_vec.x, shift_vec.y, shift_vec.z);
    // set the translation
    transformation_matrix.block<3,1>(0,3) = -shift_eigen;

    return transformation_matrix;
}

VertexData<Vector3> apply_trans_to_positions(
    const VertexData<Vector3>& positions, const Eigen::Matrix4d trans_mat
){
    VertexData<Vector3> new_positions(*positions.getMesh());
    for (Vertex v: positions.getMesh()->vertices()){
        Eigen::Vector4d pos_eigen(positions[v].x, positions[v].y, positions[v].z, 1.0);
        Eigen::Vector4d new_pos_eigen = trans_mat * pos_eigen;
        new_positions[v] = Vector3({new_pos_eigen[0], new_pos_eigen[1], new_pos_eigen[2]});
    }
    return new_positions;
}
