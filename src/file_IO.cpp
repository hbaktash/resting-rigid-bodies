#include "file_IO.h"


void validate_path(fs::path &p){
    // check if path has a parent directory that exists
    if (p.has_parent_path()){
        if (!fs::exists(p.parent_path())) {
            std::cerr << "Creating directory: " << p.parent_path() << "\n";
            fs::create_directories(p.parent_path());
        }
    }
}

nlohmann::json trans_mat_to_json(const Eigen::Matrix4d& M) {
    nlohmann::json mat = nlohmann::json::array();
    for (int r = 0; r < 4; ++r) {
        nlohmann::json row = nlohmann::json::array();
        for (int c = 0; c < 4; ++c) row.push_back(M(r,c));
        mat.push_back(row);
    }
    return mat;
}

bool save_trans_mats_and_orientations_to_file(
    const std::vector<Eigen::Matrix4d>& transformation_matrices,
    const std::vector<Vector3>& orientations, 
    const std::string& filename
) {
    fs::path out_path(filename);
    validate_path(out_path);

    nlohmann::json j;
    if (!orientations.empty()) {
        j["initial_orientation"] = {orientations[0].x, orientations[0].y, orientations[0].z};
        j["final_orientation"] = {orientations.back().x, orientations.back().y, orientations.back().z};
    }
    j["step_count"] = transformation_matrices.size();
    nlohmann::json steps = nlohmann::json::array();
    for (size_t i = 0; i < transformation_matrices.size(); ++i) {
        const auto& M = transformation_matrices[i];
        // 4x4 matrix -> nested array
        nlohmann::json mat = trans_mat_to_json(M);
        const auto& o = orientations[i];
        steps.push_back({
            {"index", (int)i},
            {"matrix", mat},
            {"orientation", {o.x, o.y, o.z}}
        });
    }
    j["steps"] = steps;

    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Error: Could not open output file: " << filename << "\n";
        return false;
    }
    ofs << j.dump(2) << "\n";
    std::cout << "Drop JSON written to " << filename << "\n";
    return true;
}

bool save_probabilities_to_file(
    Eigen::VectorXd face_region_areas, 
    Eigen::VectorX<Vector3> face_normals, 
    std::vector<Eigen::Matrix4d> stable_trans_mats,
    const std::string& filename
) {
    if (face_region_areas.size() != face_normals.size()) {
        throw std::logic_error("face_region_areas and face_normals must have the same size.");
    }
    fs::path out_path(filename);
    validate_path(out_path);

    nlohmann::json j;
    int stable_cnt = 0;
    nlohmann::json faces = nlohmann::json::array();
    for (size_t i = 0; i < face_region_areas.size(); ++i) {
        double regionArea = face_region_areas[i];
        if (regionArea > 0) {
            Vector3 n = face_normals[i];
            double prob = regionArea / (4.0 * M_PI);
            faces.push_back({
                {"hull_face_index", (int)i},
                {"face_normal", {n.x, n.y, n.z}},
                {"transformation_matrix", trans_mat_to_json(stable_trans_mats[stable_cnt++])},
                {"probability", prob}
            });
        }
    }
    j["stable_faces"] = faces;
    j["stable_face_count"] = faces.size();

    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Error: Could not open output file: " << filename << "\n";
        return false;
    }
    ofs << j.dump(2) << "\n";
    std::cout << "Probabilities JSON written to " << filename << "\n";
    return true;
}

bool save_convex_dice_params(const ConvexDiceParams& params, const std::string& filename) {
    fs::path out_path(filename);
    validate_path(out_path);
    
    nlohmann::json j;
    
    // Optimization parameters
    j["optimization"]["bary_reg"] = params.bary_reg;
    j["optimization"]["coplanar_reg"] = params.coplanar_reg;
    j["optimization"]["cluster_distance_reg"] = params.cluster_distance_reg;
    j["optimization"]["unstable_attraction_thresh"] = params.unstable_attraction_thresh;
    j["optimization"]["fair_sides_count"] = params.fair_sides_count;
    j["optimization"]["dice_energy_step"] = params.dice_energy_step;
    j["optimization"]["dice_search_decay"] = params.dice_search_decay;
    j["optimization"]["DE_step_count"] = params.DE_step_count;
    j["optimization"]["frozen_G"] = params.frozen_G;
    j["optimization"]["use_autodiff_for_dice_grad"] = params.use_autodiff_for_dice_grad;
    j["optimization"]["update_with_max_prob_face"] = params.update_with_max_prob_face;
    j["optimization"]["adaptive_reg"] = params.adaptive_reg;
    
    // Sobolev parameters
    j["sobolev"]["do_sobolev_dice_grads"] = params.do_sobolev_dice_grads;
    j["sobolev"]["sobolev_lambda"] = params.sobolev_lambda;
    j["sobolev"]["sobolev_lambda_decay"] = params.sobolev_lambda_decay;
    j["sobolev"]["sobolev_p"] = params.sobolev_p;
    
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Error: Could not open parameter file for writing: " << filename << "\n";
        return false;
    }
    
    ofs << j.dump(2) << "\n";
    std::cout << "Parameters saved to: " << filename << "\n";
    return true;
}

bool load_convex_dice_params(ConvexDiceParams& params, const std::string& filename) {
    try {
        std::ifstream ifs(filename);
        if (!ifs.is_open()) {
            std::cerr << "Error: Cannot open parameter file: " << filename << std::endl;
            return false;
        }
        
        nlohmann::json j;
        ifs >> j;
        ifs.close();
        
        // Load optimization parameters
        if (j.contains("optimization")) {
            auto opt = j["optimization"];
            if (opt.contains("bary_reg")) params.bary_reg = opt["bary_reg"];
            if (opt.contains("coplanar_reg")) params.coplanar_reg = opt["coplanar_reg"];
            if (opt.contains("cluster_distance_reg")) params.cluster_distance_reg = opt["cluster_distance_reg"];
            if (opt.contains("unstable_attraction_thresh")) params.unstable_attraction_thresh = opt["unstable_attraction_thresh"];
            if (opt.contains("fair_sides_count")) params.fair_sides_count = opt["fair_sides_count"];
            if (opt.contains("dice_energy_step")) params.dice_energy_step = opt["dice_energy_step"];
            if (opt.contains("dice_search_decay")) params.dice_search_decay = opt["dice_search_decay"];
            if (opt.contains("DE_step_count")) params.DE_step_count = opt["DE_step_count"];
            if (opt.contains("frozen_G")) params.frozen_G = opt["frozen_G"];
            if (opt.contains("use_autodiff_for_dice_grad")) params.use_autodiff_for_dice_grad = opt["use_autodiff_for_dice_grad"];
            if (opt.contains("update_with_max_prob_face")) params.update_with_max_prob_face = opt["update_with_max_prob_face"];
            if (opt.contains("adaptive_reg")) params.adaptive_reg = opt["adaptive_reg"];
        }
        
        // Load Sobolev parameters
        if (j.contains("sobolev")) {
            auto sob = j["sobolev"];
            if (sob.contains("do_sobolev_dice_grads")) params.do_sobolev_dice_grads = sob["do_sobolev_dice_grads"];
            if (sob.contains("sobolev_lambda")) params.sobolev_lambda = sob["sobolev_lambda"];
            if (sob.contains("sobolev_lambda_decay")) params.sobolev_lambda_decay = sob["sobolev_lambda_decay"];
            if (sob.contains("sobolev_p")) params.sobolev_p = sob["sobolev_p"];
        }
        
        std::cout << "Parameters loaded from: " << filename << std::endl;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error loading parameters: " << e.what() << std::endl;
        return false;
    }
}