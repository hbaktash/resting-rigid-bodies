#pragma once

#include <filesystem>
#include <fstream> 
#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <vector>
#include <string>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace fs = std::filesystem;
using Vector3 = geometrycentral::Vector3;


void validate_path(fs::path &p);

nlohmann::json trans_mat_to_json(const Eigen::Matrix4d& M);

bool save_trans_mats_and_orientations_to_file(
    const std::vector<Eigen::Matrix4d>& transformation_matrices,
    const std::vector<Vector3>& orientations, 
    const std::string& filename
);


bool save_probabilities_to_file(
    Eigen::VectorXd face_region_areas, 
    Eigen::VectorX<Vector3> face_normals, 
    std::vector<Eigen::Matrix4d> stable_trans_mats,
    const std::string& filename
);

// Parameter file I/O for convex dice optimization
struct ConvexDiceParams {
    // Optimization parameters
    float bary_reg = 0.1f;
    float coplanar_reg = 0.0f;
    float cluster_distance_reg = 0.0f;
    float unstable_attraction_thresh = 0.1f;
    int fair_sides_count = 6;
    float dice_energy_step = 0.01f;
    float dice_search_decay = 0.98f;
    int DE_step_count = 40;
    bool frozen_G = false;
    bool use_autodiff_for_dice_grad = true;
    bool update_with_max_prob_face = true;
    bool adaptive_reg = false;
    
    // Sobolev parameters
    bool do_sobolev_dice_grads = false;
    float sobolev_lambda = 2.0f;
    float sobolev_lambda_decay = 0.8f;
    int sobolev_p = 2;
};

bool save_convex_dice_params(const ConvexDiceParams& params, const std::string& filename);
bool load_convex_dice_params(ConvexDiceParams& params, const std::string& filename);
