#pragma once

#include "boundary_tools.h"
#include "utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


std::vector<std::pair<Vector3, double>> 
normal_prob_assignment(std::string shape_name);

std::vector<std::pair<Vector3, double>> 
normal_prob_assignment_fair(Forward3DSolver *tmp_solver, size_t dice_side_count);

FaceData<double> 
manual_stable_only_face_prob_assignment(Forward3DSolver *tmp_solver, std::vector<std::pair<Vector3, double>> normal_prob_pairs);

std::vector<std::tuple<std::vector<Face>, double, Vector3>> 
manual_clustered_face_prob_assignment(Forward3DSolver *tmp_solver, std::vector<std::pair<Vector3, double>> normal_prob_pairs);


// reassign cluster face normals
std::vector<std::pair<Vector3, double>> 
update_normal_prob_assignment(Forward3DSolver *tmp_solver,
                              std::vector<std::tuple<std::vector<Face>, double, Vector3>> clustered_face_normals,
                              bool take_max_prob_face);

