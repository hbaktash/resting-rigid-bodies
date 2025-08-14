#include "unsupported/Eigen/EulerAngles"

#include "forward3D.h"
#include "geometry_utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


std::pair<std::vector<Eigen::Matrix4d>, std::vector<Vector3>> 
generate_transformations_for_orientation_sequence(Vector3 initial_orientation,
    Forward3DSolver* forwardSolver, Vector3 floor_vec, double goal_angle_step);


std::pair<std::vector<Eigen::Matrix4d> , std::vector<Vector3>>
QS_trail_to_global_transformations(
    Forward3DSolver* forwardSolver, Vector3 initial_orientation, 
    Vector3 down_vec, double min_angle_step
);


VertexData<Vector3> orient_mesh_to_unit_vec(
    Vector3 orientation, Vector3 down_vec,
    ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry);

// get trans mat for orientating mesh to unit vector
Eigen::Matrix4d trans_mat_for_orientation(
    Vector3 orientation, Vector3 down_vec,
    ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry);

VertexData<Vector3> apply_trans_to_positions(
    const VertexData<Vector3>& positions, const Eigen::Matrix4d trans_mat
);
