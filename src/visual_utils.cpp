/************************************************************************
*
* ADOBE CONFIDENTIAL
* ___________________
*
* Copyright [first year code created] Adobe
* All Rights Reserved.
*
* NOTICE: All information contained herein is, and remains
* the property of Adobe and its suppliers, if any. The intellectual
* and technical concepts contained herein are proprietary to Adobe
* and its suppliers and are protected by all applicable intellectual
* property laws, including trade secret and copyright laws.
* Dissemination of this information or reproduction of this material
* is strictly forbidden unless prior written permission is obtained
* from Adobe.
*************************************************************************
*/
#include "visual_utils.h"


void draw_arc_on_sphere(Vector3 p1, Vector3 p2, Vector3 center, double radius, size_t seg_count, size_t edge_ind, polyscope::SurfaceMesh* hosting_psMesh, 
                        double radi_scale, glm::vec3 color, 
                        float arc_curve_radi, 
                        glm::vec3 default_arc_color,
                        glm::vec3 default_patch_arc_color,
                        glm::vec3 snail_trail_color){
// p1, p2 just represent normal vectors
  if (norm(p1) > 1.01 || norm(p2) > 1.01)
    polyscope::warning("wtf? p1, p2 norm larger than 1");
  
  std::vector<std::array<size_t, 2>> edgeInds;
  std::vector<Vector3> positions;
  double sqrt_radi = sqrt(radius);
  // walk on p1-p2 segment
  Vector3 curr_point = p1,
          forward_vec = (p2-p1)/(double)seg_count;
  Vector3 next_point = curr_point + forward_vec;
  Vector3 curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center ,
          next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
  positions.push_back(curr_point_on_sphere);
  for (size_t i = 0; i < seg_count; i++){
    // add to positions list
    curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center ,
    next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
    positions.push_back(next_point_on_sphere);
    // add segment indices
    edgeInds.push_back({i, i+1});

    // update points
    curr_point = next_point;
    next_point += forward_vec;
  }
  polyscope::SurfaceGraphQuantity* psArcCurve = hosting_psMesh->addSurfaceGraphQuantity("Arc curve " + std::to_string(edge_ind), positions, edgeInds);
  psArcCurve->setRadius(arc_curve_radi * radi_scale, false);
  if (color.x == -1.){
    if (edge_ind < 100)
      psArcCurve->setColor(default_arc_color);
    else if (edge_ind < 200)
      psArcCurve->setColor(default_patch_arc_color);
    else 
      psArcCurve->setColor(snail_trail_color);
  }
  else {
    psArcCurve->setColor(color);
  }
  psArcCurve->setEnabled(true);
}



void draw_arc_network_on_sphere(std::vector<std::pair<size_t, size_t>> edge_inds_,
                                std::vector<Vector3> positions_,
                                Vector3 center, double radius, size_t seg_count, std::string title, 
                                polyscope::SurfaceMesh* hosting_psMesh, 
                                double radi_scale, glm::vec3 color, 
                                float arc_curve_radi){

  std::vector<std::array<size_t, 2>> edgeInds;
  std::vector<Vector3> positions;

  double sqrt_radi = sqrt(radius);

  size_t tmp_ind = 0;
  for (std::pair<size_t, size_t> pair: edge_inds_){
    size_t i1 = pair.first, 
           i2 = pair.second;
    Vector3 p1 = positions_[i1],
            p2 = positions_[i2];
    // walk on p1-p2 segment
    Vector3 curr_point = p1,
            forward_vec = (p2-p1)/(double)seg_count;
    Vector3 next_point = curr_point + forward_vec;
    Vector3 curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center ,
            next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
    positions.push_back(curr_point_on_sphere);
    for (size_t i = 0; i < seg_count; i++){
      // add to positions list
      curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center ,
      next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
      positions.push_back(next_point_on_sphere);
      // add segment indices
      edgeInds.push_back({tmp_ind, tmp_ind+1});
      tmp_ind++;

      // update points
      curr_point = next_point;
      next_point += forward_vec;
    }
    tmp_ind++;
    // break;
  }
  printf(" in arc net: poses %d, edges %d\n", positions.size(), edgeInds.size());
  polyscope::SurfaceGraphQuantity* psArcCurve = hosting_psMesh->addSurfaceGraphQuantity("Arc curves " + title, positions, edgeInds);
  psArcCurve->setRadius(arc_curve_radi * radi_scale, false);
  psArcCurve->setColor(color);
  psArcCurve->setEnabled(true);
}