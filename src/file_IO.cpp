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