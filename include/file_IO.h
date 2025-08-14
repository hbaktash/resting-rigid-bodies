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
