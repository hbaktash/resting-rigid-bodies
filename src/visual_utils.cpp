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


// static stuff

// arc stuff
void draw_arc_on_sphere(Vector3 p1, Vector3 p2, Vector3 center, double radius, 
                        size_t seg_count, size_t edge_ind, 
                        double radi_scale, glm::vec3 color, 
                        float arc_curve_radi){
// p1, p2 just represent normal vectors
  if (norm(p1) > 1.01)
    polyscope::warning("wtf? p1 norm larger than 1");
  if (norm(p2) > 1.01)
    polyscope::warning("wtf? p2 norm larger than 1");

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
  auto psArcCurveNet = polyscope::registerCurveNetwork("Arc curve " + std::to_string(edge_ind), positions, edgeInds);
  psArcCurveNet->setRadius(arc_curve_radi * radi_scale, false);
  psArcCurveNet->setColor(color);
  psArcCurveNet->setEnabled(true);
}



void draw_arc_network_on_sphere(std::vector<std::pair<size_t, size_t>> edge_inds_,
                                std::vector<Vector3> positions_,
                                Vector3 center, double radius, size_t seg_count, std::string title, 
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
    seg_count = 3 + (int)((p1 - p2).norm()*9);
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
  // printf(" in arc net: poses %d, edges %d\n", positions.size(), edgeInds.size());
  auto psArcCurveNet = polyscope::registerCurveNetwork("Arc curves " + title, positions, edgeInds);
  psArcCurveNet->setRadius(arc_curve_radi * radi_scale, false);
  psArcCurveNet->setColor(color);
  psArcCurveNet->setEnabled(true);
}


void draw_arc_network_on_lifted_suface(std::vector<std::pair<size_t, size_t>> edge_inds_,
                                std::vector<Vector3> positions_,
                                Forward3DSolver &forward_solver,
                                Vector3 center, double radius, size_t seg_count, std::string title, 
                                glm::vec3 color, 
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
    Vector3 curr_point_on_sphere = normalize(curr_point),// * sqrt_radi,// + center ,
            next_point_on_sphere = normalize(next_point);// * sqrt_radi;// + center;
    Vector3 curr_point_on_surface = 
        center + 
        curr_point_on_sphere * sqrt(radius + forward_solver.height_function(curr_point_on_sphere));
    positions.push_back(curr_point_on_surface);
    for (size_t i = 0; i < seg_count; i++){
      // add to positions list
      curr_point_on_sphere = normalize(curr_point),
      next_point_on_sphere = normalize(next_point);
      Vector3 next_point_on_surface = 
        center + 
        next_point_on_sphere * sqrt(radius + forward_solver.height_function(next_point_on_sphere));
      positions.push_back(next_point_on_surface);
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
  // printf(" in arc net: poses %d, edges %d\n", positions.size(), edgeInds.size());
  auto psArcCurveNet = polyscope::registerCurveNetwork("Arc curves " + title, positions, edgeInds);
  psArcCurveNet->setRadius(arc_curve_radi, false);
  psArcCurveNet->setColor(color);
  psArcCurveNet->setEnabled(true);  
}


std::pair<std::vector<std::pair<size_t, size_t>>,std::vector<Vector3>> 
build_and_draw_stable_patches_on_gauss_map(BoundaryBuilder* boundary_builder, 
                                          Vector3 center, double radius, size_t seg_count,
                                          bool on_height_surface){
  std::vector<Vector3> boundary_normals(BoundaryNormal::counter);
  std::set<std::pair<size_t, size_t>> drawn_pairs;
  // printf("showing boundary patches\n");
  // printf("  building pairs\n ");
  for (Edge e: boundary_builder->forward_solver->hullMesh->edges()){
    for (BoundaryNormal *tmp_bnd_normal: boundary_builder->edge_boundary_normals[e]){
      if (tmp_bnd_normal != nullptr){
        Vector3 tmp_normal = tmp_bnd_normal->normal;
        boundary_normals[tmp_bnd_normal->index] = tmp_normal; // + gm_shift;
        for (BoundaryNormal *neigh_bnd_normal: tmp_bnd_normal->neighbors){
          Vector3 neigh_normal = neigh_bnd_normal->normal;
          boundary_normals[neigh_bnd_normal->index] = neigh_normal;
          std::pair<size_t, size_t> tmp_pair = {tmp_bnd_normal->index, neigh_bnd_normal->index};
          if (drawn_pairs.find(tmp_pair) == drawn_pairs.end()){
            drawn_pairs.insert(tmp_pair);
            // printf("here!!! %d, %d \n", tmp_pair.first, tmp_pair.second);
          }
        }
      }
    }
  }
  glm::vec3 arc_color = glm::vec3({1.,0,0}); // default color
  std::vector<std::pair<size_t, size_t>> ind_pairs_vector;
  // Using vector::assign
  ind_pairs_vector.assign(drawn_pairs.begin(), drawn_pairs.end());
  // printf("  drawing the arc network on GM\n ");
  draw_arc_network_on_sphere(ind_pairs_vector, boundary_normals, 
                            center, radius, seg_count, 
                            "region boundaries", 1., arc_color);
  if (on_height_surface){
    // printf("  drawing the arc network on height surface\n ");
    draw_arc_network_on_lifted_suface(ind_pairs_vector, boundary_normals, *boundary_builder->forward_solver, 
                                      center, 0., seg_count, 
                                      "region boundaries ", arc_color);
  }
  return {ind_pairs_vector, boundary_normals};
}


void VisualUtils::draw_edge_arcs_on_gauss_map(){
  //    add arc per edge
  std::vector<std::pair<size_t, size_t>> non_singular_edge_inds, stable_only_edge_inds ,both_edge_inds, all_edge_inds;
  size_t nFaces = forwardSolver->hullMesh->nFaces();
  std::vector<Vector3> positions(nFaces);
  for (Edge e: forwardSolver->hullMesh->edges()){
    Face f1 = e.halfedge().face(),
         f2 = e.halfedge().twin().face();
    Vector3 n1 = forwardSolver->hullGeometry->faceNormal(f1),
            n2 = forwardSolver->hullGeometry->faceNormal(f2);
    positions[f1.getIndex()] = n1;
    positions[f2.getIndex()] = n2;
    // draw with polyscope
    if (color_arcs) {
      glm::vec3 arc_color = glm::vec3({-1.,0,0}); // default color
      if (forwardSolver->edge_is_stable(e) && forwardSolver->edge_is_stablizable(e)){
        both_edge_inds.push_back({f1.getIndex(), f2.getIndex()});
        arc_color = both_edge_color;
      }
      else if (forwardSolver->edge_is_stable(e)){
        stable_only_edge_inds.push_back({f1.getIndex(), f2.getIndex()});
        arc_color = stable_edge_color;
      }
      else if (draw_unstable_edge_arcs) {
        if (forwardSolver->edge_is_stablizable(e)){
          non_singular_edge_inds.push_back({f1.getIndex(), f2.getIndex()});
          arc_color = stabilizable_edge_color;
        }
        all_edge_inds.push_back({f1.getIndex(), f2.getIndex()});
      }
    } // else; draw all arcs with black
    else {
      all_edge_inds.push_back({f1.getIndex(), f2.getIndex()});
    }
  }
  if (color_arcs){
    draw_arc_network_on_sphere(both_edge_inds, positions, center, gm_radi, arcs_seg_count,
                              "saddle edge arcs", 1., both_edge_color);
    draw_arc_network_on_sphere(stable_only_edge_inds, positions, center, gm_radi, arcs_seg_count,
                              "singular edge arcs", 1., stable_edge_color);
    draw_arc_network_on_sphere(all_edge_inds, positions, center, gm_radi, arcs_seg_count,
                              "non-singular edge arcs", 1., stabilizable_edge_color);
  }
  else {
    draw_arc_network_on_sphere(all_edge_inds, positions, center, gm_radi, arcs_seg_count,
                              "all edge arcs", 1., stabilizable_edge_color);
  }
}


void VisualUtils::draw_stable_vertices_on_gauss_map(){
  std::vector<Vector3> stable_vertices, hidden_stable_vertices;
  for (Vertex v: forwardSolver->hullMesh->vertices()){
    if (forwardSolver->vertex_is_stablizable(v)){
      stable_vertices.push_back(normalize(forwardSolver->hullGeometry->inputVertexPositions[v] - forwardSolver->get_G())+center);
    }
    else {
      hidden_stable_vertices.push_back(normalize(forwardSolver->hullGeometry->inputVertexPositions[v] - forwardSolver->get_G())+center);
    }
  }
  polyscope::PointCloud* stable_vertices_gm_pc = polyscope::registerPointCloud("stable Vertices Normals", 
                                                                               stable_vertices);
  stable_vertices_gm_pc->setPointRadius(gm_pt_radi, false);
  stable_vertices_gm_pc->setPointColor({0.1,0.9,0.1});
  stable_vertices_gm_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

  // hidden stable vertices
  if (show_hidden_stable_vertex_normals){
    polyscope::PointCloud* hidden_stable_vertices_gm_pc = polyscope::registerPointCloud(
                                                          "hidden stable Vertices Normals", 
                                                          hidden_stable_vertices);
    hidden_stable_vertices_gm_pc->setPointRadius(gm_pt_radi, false);
    hidden_stable_vertices_gm_pc->setPointColor({0.03,0.4,0.03});
    hidden_stable_vertices_gm_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  }
  else {
    if (polyscope::hasPointCloud("hidden stable Vertices Normals"))
      polyscope::removePointCloud("hidden stable Vertices Normals");
  }
}



void VisualUtils::draw_stable_face_normals_on_gauss_map(){
  std::vector<Vector3> stable_face_normals, face_normal_points;
  for (Face f: forwardSolver->hullMesh->faces()){
    Vector3 normal_pos_on_gm = forwardSolver->hullGeometry->faceNormal(f) + center;
    face_normal_points.push_back(normal_pos_on_gm);
    if (forwardSolver->face_is_stable(f)){
      stable_face_normals.push_back(normal_pos_on_gm);
    }
  }
  polyscope::PointCloud* stable_face_normals_pc = polyscope::registerPointCloud("stable Face Normals", stable_face_normals);
  stable_face_normals_pc->setPointRadius(gm_pt_radi*1.1, false);
  stable_face_normals_pc->setPointColor({0.9,0.1,0.1});
  stable_face_normals_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere); 

  polyscope::PointCloud* face_normals_pc = polyscope::registerPointCloud("Face Normals", face_normal_points);
  face_normals_pc->setPointRadius(gm_pt_radi, false);
  face_normals_pc->setPointColor({0.9,0.9,0.9});
  face_normals_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);
}



void VisualUtils::plot_height_function(){
  VertexData<double> heigh_func(*sphere_mesh);
  VertexData<Vector3> lifted_poses(*sphere_mesh);
  // printf("here! sphre mesh sizE: %d \n", sphere_mesh->nVertices());
  for (Vertex v: sphere_mesh->vertices()){
    Vector3 pos = sphere_geometry->inputVertexPositions[v];
    double height = forwardSolver->height_function(pos.normalize());
    heigh_func[v] = height;
    lifted_poses[v] = sqrt(height) * pos.normalize() + center; // 1. + height
  }
  polyscope::SurfaceVertexScalarQuantity *height_func_vis = 
              gm_sphere_mesh->addVertexScalarQuantity(" height function", heigh_func);
  // height_func_vis->setEnabled(true);

  polyscope::SurfaceMesh* height_function_mesh = polyscope::registerSurfaceMesh("height surface_func", 
                                                  lifted_poses, 
                                                  sphere_mesh->getFaceVertexList());
  height_function_mesh->setSmoothShade(true);
  polyscope::SurfaceVertexScalarQuantity *height_func_vis_surface = 
              height_function_mesh->addVertexScalarQuantity(" height function surface", heigh_func);
  // height_func_vis_surface->setEnabled(true);
}


void VisualUtils::draw_gauss_map(){
  // just draw the sphere next to the main surface
  // std::vector<Vector3> sphere_pos = {gm_shift};
  // gauss_map_pc = polyscope::registerPointCloud("Gauss Map", sphere_pos);
  // gauss_map_pc->setPointColor({0.74,0.7,0.9});
  // gauss_map_pc->setPointRadius(gm_radi, false);
  // gauss_map_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  VertexData<Vector3> shifted_poses(*sphere_mesh);
  for (Vertex v: sphere_mesh->vertices()){
    Vector3 pos = sphere_geometry->inputVertexPositions[v];
    shifted_poses[v] = pos.normalize() + center;
  }
  gm_sphere_mesh = polyscope::registerSurfaceMesh("gm_sphere_mesh", 
                                                  shifted_poses, 
                                                  sphere_mesh->getFaceVertexList());
  gm_sphere_mesh->setSmoothShade(true);
  gm_sphere_mesh->setSurfaceColor({0.74,0.7,0.9});
  // surface and colormap
  plot_height_function();
  
  // face normals on Gauss map
  draw_stable_face_normals_on_gauss_map();

  // point cloud for stable vertices
  draw_stable_vertices_on_gauss_map();

  // arcs for edge-normals set
  draw_edge_arcs_on_gauss_map();

  // edge equilibria
  show_edge_equilibria_on_gauss_map();
}


void VisualUtils::show_edge_equilibria_on_gauss_map(){
  std::vector<Vector3> edge_equilibria_points, stabilizable_edge_equilibria_points, stable_edge_equilibria_points;
  Vector3 G = forwardSolver->get_G();
  for (Edge e: forwardSolver->hullMesh->edges()){
    if (forwardSolver->edge_is_stablizable(e) && forwardSolver->edge_is_stable(e)){
      Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
      Vector3 A = forwardSolver->hullGeometry->inputVertexPositions[v1], B = forwardSolver->hullGeometry->inputVertexPositions[v2];
      Vector3 GB = B - G,
              AB = B - A;
      Vector3 ortho_g = GB - AB*dot(AB, GB)/dot(AB,AB);
      edge_equilibria_points.push_back(normalize(ortho_g) + center);
    }
    else if (forwardSolver->edge_is_stablizable(e)){
      Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
      Vector3 A = forwardSolver->hullGeometry->inputVertexPositions[v1], B = forwardSolver->hullGeometry->inputVertexPositions[v2];
      Vector3 GB = B - G,
              AB = B - A;
      Vector3 ortho_g = GB - AB*dot(AB, GB)/dot(AB,AB);
      stabilizable_edge_equilibria_points.push_back(normalize(ortho_g) + center);
    }
    else if (forwardSolver->edge_is_stable(e)){
      Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
      Vector3 A = forwardSolver->hullGeometry->inputVertexPositions[v1], B = forwardSolver->hullGeometry->inputVertexPositions[v2];
      Vector3 GB = B - G,
              AB = B - A;
      Vector3 ortho_g = GB - AB*dot(AB, GB)/dot(AB,AB);
      stable_edge_equilibria_points.push_back(normalize(ortho_g) + center);
    }
  }
  polyscope::PointCloud* edge_equilibria_pc = polyscope::registerPointCloud("Edge equilibria", edge_equilibria_points);
  edge_equilibria_pc->setPointRadius(gm_pt_radi, false);
  edge_equilibria_pc->setPointColor({0.2, 0.2, 0.9});
  edge_equilibria_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

  if (draw_stable_g_vec_for_unstable_edge_arcs){
    polyscope::PointCloud* stabilizable_edge_pc = polyscope::registerPointCloud(
                      "Stabilizable Edge equilibria", stabilizable_edge_equilibria_points);
    stabilizable_edge_pc->setPointRadius(gm_pt_radi, false);
    stabilizable_edge_pc->setPointColor(stabilizable_edge_color);
    stabilizable_edge_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

    polyscope::PointCloud* stable_edge_pc = polyscope::registerPointCloud("stable Edge equilibria", stable_edge_equilibria_points);
    stable_edge_pc->setPointRadius(gm_pt_radi, false);
    stable_edge_pc->setPointColor(stable_edge_color);
    stable_edge_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  }
  else{
    if (polyscope::hasPointCloud("stable Edge equilibria"))
      polyscope::removePointCloud("stable Edge equilibria");
    if (polyscope::hasPointCloud("Stabilizable Edge equilibria"))
      polyscope::removePointCloud("Stabilizable Edge equilibria");
  }
}

void VisualUtils::draw_guess_pc(std::vector<std::pair<size_t, size_t>> neigh_inds, 
                                std::vector<Vector3> boundary_normals){
  std::vector<Vector3> positions;
  std::vector<std::array<size_t, 2>> edgeInds;
  for (Vector3 bnd_normal: boundary_normals){
    Vector3 pos_on_mesh;
    double tmp_t = -1;
    for (Face f: forwardSolver->hullMesh->faces()){
      std::vector<Vector3> face_poses;
      for (Vertex v: f.adjacentVertices())
        face_poses.push_back(forwardSolver->hullGeometry->inputVertexPositions[v]);
      tmp_t = ray_intersect(forwardSolver->get_G(), bnd_normal, face_poses);
      if (tmp_t != -1)
        break;
    }
    // on surface
    positions.push_back(forwardSolver->get_G() + tmp_t * bnd_normal);
    // easy case
    // positions.push_back(bnd_normal + shift);
  }
  polyscope::PointCloud* test_pc = polyscope::registerPointCloud("test point cloud", positions);
  test_pc->setPointColor({1.,0.,0.});
  // test_pc->setEnabled(true);
  test_pc->setPointRadius(gm_pt_radi * 0.7, false);
  test_pc->setEnabled(true);

  std::vector<std::array<size_t, 2>> edge_inds;
  for (std::pair<size_t, size_t> p1: neigh_inds)
    edge_inds.push_back({p1.first, p1.second});
  auto stable_graph_PSnet = 
        polyscope::registerCurveNetwork("stable regions on polyhedra", positions, edge_inds);
  stable_graph_PSnet->setRadius(stable_edge_radi/2., true);
  stable_graph_PSnet->setColor(stable_edge_color);
  stable_graph_PSnet->setEnabled(true);
  
}


// visualize center of mass
void VisualUtils::draw_G() {
  std::vector<Vector3> G_position = {forwardSolver->get_G()};
  psG = polyscope::registerPointCloud("Center of Mass", G_position);
  psG->setPointColor({0., 0., 0.});
  psG->setPointRadius(G_radi, false);
}



// polyhedra domain

// void visualize_edge_probabilities(){}
void VisualUtils::visualize_stable_vertices(){
  std::vector<Vector3> positions;// = {forwardSolver->hullGeometry->inputVertexPositions[v]};
  for (Vertex v: forwardSolver->hullMesh->vertices()){
    if (forwardSolver->vertex_is_stablizable(v)){
      positions.push_back(forwardSolver->hullGeometry->inputVertexPositions[v]);
    }
  }
  polyscope::PointCloud* psCloud = polyscope::registerPointCloud("stabilizable vertices", positions);
  // set some options
  psCloud->setPointColor({0.1, .9, .1});
  psCloud->setPointRadius(gm_pt_radi);
  psCloud->setPointRenderMode(polyscope::PointRenderMode::Sphere);
}




void VisualUtils::visualize_edge_stability(){
  std::vector<std::array<size_t, 2>> stable_edgeInds, stablilizable_edgeInds, both_edgeInds;
  std::vector<Vector3> stable_positions, stablilizable_positions, both_positions;
  size_t stable_counter = 0, stablizable_counter = 0, both_counter = 0;
  for (Edge e: forwardSolver->hullMesh->edges()){
    Vector3 p1 = forwardSolver->hullGeometry->inputVertexPositions[e.firstVertex()],
            p2 = forwardSolver->hullGeometry->inputVertexPositions[e.secondVertex()];
    size_t flag = 0;
    if (forwardSolver->edge_is_stable(e)){
      stable_positions.push_back(p1); stable_positions.push_back(p2);
      stable_edgeInds.push_back({stable_counter, stable_counter + 1});
      stable_counter += 2;
      flag++;
    }
    if (forwardSolver->edge_is_stablizable(e)){
      stablilizable_positions.push_back(p1); stablilizable_positions.push_back(p2);
      stablilizable_edgeInds.push_back({stablizable_counter, stablizable_counter + 1});
      stablizable_counter += 2;
      flag++;
    }
    if (flag == 2){
      both_positions.push_back(p1); both_positions.push_back(p2);
      both_edgeInds.push_back({both_counter, both_counter + 1});
      both_counter += 2;
    }
  }
  auto psStableEdgesNet =  polyscope::registerCurveNetwork("stable Edges", stable_positions, stable_edgeInds);
  psStableEdgesNet->setRadius(stable_edge_radi, true);
  psStableEdgesNet->setColor(stable_edge_color);
  psStableEdgesNet->setEnabled(true);
  auto psStablilizableEdgesNet =  polyscope::registerCurveNetwork("stablizable Edges", stablilizable_positions, stablilizable_edgeInds);
  psStablilizableEdgesNet->setRadius(stablizable_edge_radi, true);
  psStablilizableEdgesNet->setColor(stabilizable_edge_color);
  psStablilizableEdgesNet->setEnabled(true);
  auto psBothEdgesNet =  polyscope::registerCurveNetwork("stable && stablizable Edges", both_positions, both_edgeInds);
  psBothEdgesNet->setRadius(both_edge_radi, true);
  psBothEdgesNet->setColor(both_edge_color);
  psBothEdgesNet->setEnabled(true);
}


void VisualUtils::visualize_face_stability(){
  std::vector<std::array<double, 3>> fColor(forwardSolver->hullMesh->nFaces());
  for (Face f: forwardSolver->hullMesh->faces()){
    if (forwardSolver->face_is_stable(f))
      fColor[f.getIndex()] = {1., 0.1, 0.1};
    else
      fColor[f.getIndex()] = {0.9, 0.9, 1.};
  }
  polyscope::SurfaceMesh* psInputMesh = polyscope::getSurfaceMesh("hull mesh");
  polyscope::SurfaceFaceColorQuantity *faceQnty =  psInputMesh->addFaceColorQuantity("face stability", fColor);
  faceQnty->setEnabled(true);
}

void VisualUtils::visualize_colored_polyhedra(FaceData<Vector3> face_colors){
  VertexData<Vector3> shifted_positions(*forwardSolver->hullMesh);
  for (Vertex v: forwardSolver->hullMesh->vertices()){
    shifted_positions[v] = forwardSolver->hullGeometry->inputVertexPositions[v] + colored_shift;
  }
  // printf("hull faces: %d color size: %d \n", forwardSolver->hullMesh->nFaces(), face_colors.size());
  polyscope::SurfaceMesh* coloredPsMesh = polyscope::registerSurfaceMesh("colored polyhedra", shifted_positions, forwardSolver->hullMesh->getFaceVertexList());
  // generate random colors and color the faces
  polyscope::SurfaceFaceColorQuantity *faceQnty = coloredPsMesh->addFaceColorQuantity("random face colors", face_colors);
  faceQnty->setEnabled(true);
  // // add colors to the original polyhedra as well?
  // polyscope::SurfaceFaceColorQuantity *faceQnty2 = psMesh->addFaceColorQuantity("random face colors2", face_colors);
  // faceQnty2->setEnabled(true);
}


void VisualUtils::visualize_all_stable_orientations(){
  forwardSolver->set_uniform_G();
  forwardSolver->initialize_pre_computes();
  BoundaryBuilder *bnd_builder = new BoundaryBuilder(forwardSolver);
  bnd_builder->build_boundary_normals();

  std::vector<std::pair<Face, double>> probs;
  for (Face f: forwardSolver->hullMesh->faces())
      if (bnd_builder->face_region_area[f] > 0)
          probs.push_back({f, bnd_builder->face_region_area[f]/(4.*PI)});
  std::sort(probs.begin(), probs.end(), [] (auto a, auto b) { return a.second > b.second; });
  double floor_z = -1.;
  Vector3 floor_vec = Vector3({0,0,floor_z});
  VertexData<Vector3> rotated_poses(*forwardSolver->inputMesh);
  size_t i = 0;
  for (auto &p: probs){
    Face f = p.first;
    double prob = p.second;
    Vector3 f_normal = forwardSolver->hullGeometry->faceNormal(f);
    Vector3 rot_axis = cross(f_normal, floor_vec);
    double rot_angle = angle(f_normal, floor_vec);
    Vector3 offset = floor_vec - forwardSolver->hullGeometry->inputVertexPositions[f.halfedge().vertex()].rotateAround(rot_axis, rot_angle);
    for (Vertex v: forwardSolver->inputMesh->vertices()){
      rotated_poses[v] = forwardSolver->inputGeometry->inputVertexPositions[v].rotateAround(rot_axis, rot_angle) + offset;
    }
    double y_shift = - (double)i;
    Vector3 vis_shift({y_shift * 2., -2., 0.});
    auto tmp_ori_mesh = polyscope::registerSurfaceMesh("ori "+std::to_string(i)+" probe: " + std::to_string(prob), rotated_poses + vis_shift, forwardSolver->inputMesh->getFaceVertexList());
    tmp_ori_mesh->setSurfaceColor({0.1, 0.4, 0.02});
    i++;
  }

}