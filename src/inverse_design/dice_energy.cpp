#include "inverse_design/dice_energy.h"



void get_dice_energy_grads(Eigen::MatrixX3d hull_positions, Eigen::Vector3d G_vec, double bary_reg, double coplanar_reg, double cluster_distance_reg, double unstable_attraction_thresh,
                           Eigen::MatrixX3d &df_dv, Eigen::Vector3d &df_dG, double &dice_energy,
                           bool frozen_G, 
                           std::string policy_general, std::vector<std::pair<Vector3, double>> normal_prob_pairs, 
                           int fair_sides){
  // Eigen::MatrixX3d hull_positions = vertex_data_to_matrix(fwd_solver.hullGeometry->inputVertexPositions); 
  Forward3DSolver tmp_solver(hull_positions, G_vec, true); // indices shouldnt be shuffled here
  auto dice_energy_lambda = [&] <typename Scalar> (const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &hull_poses_G_append_vec) -> Scalar {
    // decompose flat vector to positions and center of mass; G is the last 3 elements
    Eigen::Vector3<Scalar> G_eigen = hull_poses_G_append_vec.tail(3);
    size_t flat_n = hull_poses_G_append_vec.rows();
    Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 3> > hull_poses(hull_poses_G_append_vec.head(flat_n-3).data(), flat_n/3 - 1, 3);
    return ::dice_energy<Scalar>(hull_poses, G_eigen, tmp_solver, 
								bary_reg, coplanar_reg, cluster_distance_reg, unstable_attraction_thresh,
								policy_general, normal_prob_pairs, fair_sides, true);
  };
  Eigen::VectorXd hull_poses_vec = hull_positions.reshaped();
  Eigen::VectorXd hull_poses_and_G_vec(hull_poses_vec.size() + 3);
  hull_poses_and_G_vec << hull_poses_vec, G_vec;
  
  Eigen::VectorXd dfdU_vec;
  double dice_e;
  stan::math::gradient(dice_energy_lambda, hull_poses_and_G_vec, dice_e, dfdU_vec);
  dice_energy = dice_e;
  size_t flat_n = dfdU_vec.rows();
  df_dG = dfdU_vec.tail(3);
  Eigen::Map<Eigen::MatrixXd> dfdV(dfdU_vec.head(flat_n-3).data(), flat_n/3 - 1, 3);

  // populate df_dv by mapping to original input indices
  if (frozen_G){
    df_dv = dfdV;
  }
  else {
    std::vector<Eigen::Matrix3d> dG_dv = get_COM_grads_for_convex_uniform_shape(hull_positions);
    for (size_t i = 0; i < hull_positions.rows(); i++){
      df_dv.row(i) = dfdV.row(i) + (dG_dv[i].transpose() * df_dG).transpose();
    }
  }  
}

// hull update stuff
double hull_update_line_search(Eigen::MatrixX3d dfdv, Eigen::MatrixX3d hull_positions, Eigen::Vector3d G_vec, 
                               double bary_reg, double coplanar_reg, double cluster_distance_reg, double unstable_attaction_thresh,
                               std::string policy_general, std::vector<std::pair<Vector3, double>> normal_prob_assignment, 
                               size_t dice_side_count, 
                               double step_size, double decay, bool frozen_G, size_t max_iter, double step_tol){
  
  Forward3DSolver tmp_solver(hull_positions, G_vec, true); // assuming input is convex; will be asserted internally in the constructor
  if (!frozen_G){
    tmp_solver.set_uniform_G();
    G_vec = vec32vec(tmp_solver.get_G());
  }
  tmp_solver.initialize_pre_computes();
  
  // std::cout << "getting fair dice energy for the initial hull\n";
  double min_dice_energy = dice_energy<double>(hull_positions, G_vec, 
                                                                tmp_solver, bary_reg, coplanar_reg, cluster_distance_reg, unstable_attaction_thresh,
                                                                policy_general, normal_prob_assignment, 
                                                                dice_side_count, false);
  double s = step_size; //

  bool found_smth_optimal = false;
  double tmp_dice_energy;
  int j;
  for (j = 0; j < max_iter; j++) {
    tmp_solver = Forward3DSolver(hull_positions - s * dfdv, G_vec, false); // not necessarily convex
    if (frozen_G && !G_is_inside(*tmp_solver.hullMesh, *tmp_solver.hullGeometry, tmp_solver.get_G())){
        // printf("  - G outside! \n");
        s *= decay;
        continue;
    }
    if (!frozen_G){
        tmp_solver.set_uniform_G();
        G_vec = vec32vec(tmp_solver.get_G());
    }
    tmp_solver.initialize_pre_computes();
    // re-assign since qhull inside solver reshuffles points
    Eigen::MatrixX3d tmp_hull_positions = vertex_data_to_matrix(tmp_solver.hullGeometry->inputVertexPositions);
    Forward3DSolver tmp_solver2(tmp_hull_positions, G_vec, true);
    bool verbose = false;
    if (s < step_tol){
      verbose = true;
      printf("   ---   LS step %d -----\n", j);
    }
    tmp_dice_energy = dice_energy<double>(tmp_hull_positions, G_vec,
                                            tmp_solver2, bary_reg, coplanar_reg, cluster_distance_reg, unstable_attaction_thresh,
                                            policy_general, normal_prob_assignment, 
                                            dice_side_count, verbose);

    if (tmp_dice_energy < min_dice_energy){
        found_smth_optimal = true;
        break; //  x new is good
    }
    else if(s >= step_tol)
        s *= decay;
    else
      break;
  }
  s = found_smth_optimal ? s : 0.;
  printf("line search for dice ended at iter %d, s: %.10f, \n", j, s);
  return s;
}


