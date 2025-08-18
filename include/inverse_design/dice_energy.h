#pragma once

#include "boundary_tools.h"
#include "utils.h"
#include "set"
#include "inverse_design/prob_assignment.h"
#include <stan/math.hpp>

// Template stuff for inverse design
template <typename Scalar>
 Scalar dice_energy(Eigen::MatrixX3<Scalar> hull_positions, Eigen::Vector3<Scalar> G,
                            Forward3DSolver &tmp_solver, double bary_reg, double co_planar_reg, double cluster_distance_reg, double unstable_attaction_thresh,
                            std::string policy, std::vector<std::pair<Vector3, double>> normal_prob_assignment, 
                            size_t side_count, bool verbose);
template <typename Scalar>
 Eigen::Vector3<Scalar> point_to_segment_normal(Eigen::Vector3<Scalar> P, Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B);
template <typename Scalar>
 Eigen::Vector3<Scalar> intersect_arcs(Eigen::Vector3<Scalar> v_normal, Eigen::Vector3<Scalar> R2, Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B);
template <typename Scalar>
 Scalar triangle_patch_signed_area_on_sphere(Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B, Eigen::Vector3<Scalar> C);

// Regularizers
//  --- for ManualCluster policy --- 
template <typename Scalar>
 Scalar 
single_cluster_coplanar_e(std::vector<Face> faces, 
                            FaceData<Eigen::Vector3<Scalar>> face_normals,
                            FaceData<Face> face_last_face, double unstable_attaction_thresh, 
                            bool verbose);

template <typename Scalar>
 Scalar 
single_cluster_bary_e(std::vector<Face> faces, FaceData<Face> face_last_face,
                        Eigen::MatrixX3<Scalar> hull_positions, Eigen::Vector3<Scalar> G,
                        FaceData<Eigen::Vector3<Scalar>> face_normals,
                        FaceData<std::set<Vertex>> face_region_boundary_vertices);
template <typename Scalar>
 std::pair<Eigen::Vector3<Scalar>, Eigen::Vector3<Scalar>> 
region_and_stables_barycenters(std::vector<Face> faces, FaceData<Face> face_last_face, 
                                Eigen::MatrixX3<Scalar> hull_positions, Eigen::Vector3<Scalar> G,
                                FaceData<Eigen::Vector3<Scalar>> face_normals,
                                FaceData<std::set<Vertex>> face_region_boundary_vertices);

// --- for Fair and Manual policies ---
template <typename Scalar>
 Scalar 
single_DAG_cluster_bary_e(Face stable_face, Eigen::Vector3<Scalar> stable_face_normal,
                        Eigen::MatrixX3<Scalar> hull_positions, Eigen::Vector3<Scalar> G,
                        FaceData<std::set<Vertex>> face_region_boundary_vertices);


void get_dice_energy_grads(Eigen::MatrixX3d hull_positions, Eigen::Vector3d G_vec, double bary_reg, double coplanar_reg, double cluster_distance_reg, double unstable_attraction_thresh,
    Eigen::MatrixX3d &df_dv, Eigen::Vector3d &df_dG, double &dice_energy,
    bool frozen_G, 
    std::string policy_general, std::vector<std::pair<Vector3, double>> normal_prob_pairs, 
    int fair_sides);


double hull_update_line_search(Eigen::MatrixX3d dfdv, Eigen::MatrixX3d hull_positions, Eigen::Vector3d G_vec, 
                               double bary_reg, double coplanar_reg, double cluster_distance_reg, double unstable_attaction_thresh,
                               std::string policy_general, std::vector<std::pair<Vector3, double>> normal_prob_assignment, 
                               size_t dice_side_count, double step_size, double decay, bool frozen_G, 
                               size_t max_iter, double step_tol);

#include "inverse_design/dice_energy.impl.h"