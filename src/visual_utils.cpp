/************************************************************************
*
* ADOBE CONFIDENTIAL
* ___________________
*
* (Original notice retained; replace for publication.)
*************************************************************************
*/
#include "visual_utils.h"

// arc helpers (already free)
void draw_arc_on_sphere(Vector3 p1, Vector3 p2, Vector3 center, double radius,
                        size_t seg_count, size_t edge_ind,
                        double radi_scale, glm::vec3 color,
                        float arc_curve_radi) {
  if (norm(p1) > 1.01) polyscope::warning("wtf? p1 norm larger than 1");
  if (norm(p2) > 1.01) polyscope::warning("wtf? p2 norm larger than 1");

  std::vector<std::array<size_t, 2>> edgeInds;
  std::vector<Vector3> positions;
  double sqrt_radi = sqrt(radius);

  Vector3 curr_point = p1, forward_vec = (p2 - p1) / (double)seg_count;
  Vector3 next_point = curr_point + forward_vec;
  Vector3 curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center,
          next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
  positions.push_back(curr_point_on_sphere);
  for (size_t i = 0; i < seg_count; i++) {
    curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center,
    next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
    positions.push_back(next_point_on_sphere);
    edgeInds.push_back({i, i + 1});
    curr_point = next_point;
    next_point += forward_vec;
  }
  auto psArcCurveNet =
      polyscope::registerCurveNetwork("Arc curve " + std::to_string(edge_ind), positions, edgeInds);
  psArcCurveNet->setRadius(arc_curve_radi * radi_scale, false);
  psArcCurveNet->setColor(color);
  psArcCurveNet->setEnabled(true);
}

void draw_arc_network_on_sphere(std::vector<std::pair<size_t, size_t>> edge_inds_,
                                std::vector<Vector3> positions_, Vector3 center, double radius,
                                size_t seg_count, std::string title, double radi_scale,
                                glm::vec3 color, float arc_curve_radi) {

  std::vector<std::array<size_t, 2>> edgeInds;
  std::vector<Vector3> positions;
  double sqrt_radi = sqrt(radius);

  size_t tmp_ind = 0;
  for (auto& pair : edge_inds_) {
    size_t i1 = pair.first, i2 = pair.second;
    Vector3 p1 = positions_[i1], p2 = positions_[i2];
    int tmp_seg_cnt = 3 + (int)((p1 - p2).norm() * seg_count);
    Vector3 curr_point = p1, forward_vec = (p2 - p1) / (double)tmp_seg_cnt;
    Vector3 next_point = curr_point + forward_vec;
    Vector3 curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center,
            next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
    positions.push_back(curr_point_on_sphere);
    for (size_t i = 0; i < (size_t)tmp_seg_cnt; i++) {
      curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center,
      next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
      positions.push_back(next_point_on_sphere);
      edgeInds.push_back({tmp_ind, tmp_ind + 1});
      tmp_ind++;
      curr_point = next_point;
      next_point += forward_vec;
    }
    tmp_ind++;
  }
  auto psArcCurveNet =
      polyscope::registerCurveNetwork("Arc curves " + title, positions, edgeInds);
  psArcCurveNet->setRadius(arc_curve_radi * radi_scale, false);
  psArcCurveNet->setColor(color);
  psArcCurveNet->setEnabled(true);
}

void draw_arc_network_on_lifted_suface(std::vector<std::pair<size_t, size_t>> edge_inds_,
                                       std::vector<Vector3> positions_,
                                       Forward3DSolver& forward_solver, Vector3 center,
                                       double radius, size_t seg_count, std::string title,
                                       glm::vec3 color, float arc_curve_radi) {
  std::vector<std::array<size_t, 2>> edgeInds;
  std::vector<Vector3> positions;
  double sqrt_radi = sqrt(radius);
  size_t tmp_ind = 0;
  for (auto& pair : edge_inds_) {
    size_t i1 = pair.first, i2 = pair.second;
    Vector3 p1 = positions_[i1], p2 = positions_[i2];
    Vector3 curr_point = p1, forward_vec = (p2 - p1) / (double)seg_count;
    Vector3 next_point = curr_point + forward_vec;
    Vector3 curr_point_on_sphere = normalize(curr_point),
            next_point_on_sphere = normalize(next_point);
    Vector3 curr_point_on_surface =
        center + curr_point_on_sphere *
                     sqrt(radius + forward_solver.height_function(curr_point_on_sphere));
    positions.push_back(curr_point_on_surface);
    for (size_t i = 0; i < seg_count; i++) {
      curr_point_on_sphere = normalize(curr_point),
      next_point_on_sphere = normalize(next_point);
      Vector3 next_point_on_surface =
          center + next_point_on_sphere *
                       sqrt(radius + forward_solver.height_function(next_point_on_sphere));
      positions.push_back(next_point_on_surface);
      edgeInds.push_back({tmp_ind, tmp_ind + 1});
      tmp_ind++;
      curr_point = next_point;
      next_point += forward_vec;
    }
    tmp_ind++;
  }
  auto psArcCurveNet =
      polyscope::registerCurveNetwork("Arc curves " + title, positions, edgeInds);
  psArcCurveNet->setRadius(arc_curve_radi, false);
  psArcCurveNet->setColor(color);
  psArcCurveNet->setEnabled(true);
}

std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<Vector3>>
build_and_draw_stable_patches_on_gauss_map(BoundaryBuilder* boundary_builder, Vector3 center,
                                           double radius, size_t seg_count, bool on_height_surface,
                                           double arc_curve_radi) {
  std::vector<Vector3> boundary_normals(BoundaryNormal::counter);
  std::vector<std::pair<size_t, size_t>> ind_pairs_vector;

  auto p = boundary_builder->MS_complex_edges();
  ind_pairs_vector = p.first;
  boundary_normals = p.second;

  glm::vec3 arc_color = glm::vec3({1., 0, 0});
  draw_arc_network_on_sphere(ind_pairs_vector, boundary_normals, center, radius, seg_count,
                             "region boundaries", 1., arc_color, arc_curve_radi);
  if (on_height_surface) {
    draw_arc_network_on_lifted_suface(ind_pairs_vector, boundary_normals,
                                      *boundary_builder->forward_solver, center, 0., seg_count,
                                      "region boundaries ", arc_color, arc_curve_radi);
  }
  return {ind_pairs_vector, boundary_normals};
}

// --- Converted member functions -> free functions ----------------------------

void draw_edge_arcs_on_gauss_map(Forward3DSolver* forwardSolver) {
  std::vector<std::pair<size_t, size_t>> non_singular_edge_inds, stable_only_edge_inds,
      both_edge_inds, all_edge_inds;
  size_t nFaces = forwardSolver->hullMesh->nFaces();
  std::vector<Vector3> positions(nFaces);
  for (Edge e : forwardSolver->hullMesh->edges()) {
    Face f1 = e.halfedge().face(), f2 = e.halfedge().twin().face();
    Vector3 n1 = forwardSolver->hullGeometry->faceNormal(f1),
            n2 = forwardSolver->hullGeometry->faceNormal(f2);
    positions[f1.getIndex()] = n1;
    positions[f2.getIndex()] = n2;
    if (VisualUtils::color_arcs) {
      if (forwardSolver->edge_is_stable(e) && forwardSolver->edge_is_stablizable(e)) {
        both_edge_inds.push_back({f1.getIndex(), f2.getIndex()});
      } else if (forwardSolver->edge_is_stable(e)) {
        stable_only_edge_inds.push_back({f1.getIndex(), f2.getIndex()});
      } else if (VisualUtils::draw_unstable_edge_arcs) {
        if (forwardSolver->edge_is_stablizable(e)) {
          non_singular_edge_inds.push_back({f1.getIndex(), f2.getIndex()});
        }
        all_edge_inds.push_back({f1.getIndex(), f2.getIndex()});
      }
    } else {
      all_edge_inds.push_back({f1.getIndex(), f2.getIndex()});
    }
  }
  if (VisualUtils::color_arcs) {
    draw_arc_network_on_sphere(both_edge_inds, positions, VisualUtils::center, VisualUtils::gm_radi,
                               VisualUtils::arcs_seg_count, "saddle edge arcs", 1.,
                               VisualUtils::both_edge_color);
    draw_arc_network_on_sphere(stable_only_edge_inds, positions, VisualUtils::center,
                               VisualUtils::gm_radi, VisualUtils::arcs_seg_count,
                               "singular edge arcs", 1., VisualUtils::stable_edge_color);
    draw_arc_network_on_sphere(all_edge_inds, positions, VisualUtils::center, VisualUtils::gm_radi,
                               VisualUtils::arcs_seg_count, "non-singular edge arcs", 1.,
                               VisualUtils::stabilizable_edge_color);
  } else {
    draw_arc_network_on_sphere(all_edge_inds, positions, VisualUtils::center, VisualUtils::gm_radi,
                               VisualUtils::arcs_seg_count, "all edge arcs", 1.,
                               VisualUtils::stabilizable_edge_color,
                               VisualUtils::arc_curve_radi);
  }
}

void draw_stable_vertices_on_gauss_map(Forward3DSolver* forwardSolver) {
  std::vector<Vector3> stable_vertices, hidden_stable_vertices;
  for (Vertex v : forwardSolver->hullMesh->vertices()) {
    if (forwardSolver->vertex_is_stablizable(v)) {
      stable_vertices.push_back(normalize(
                                    forwardSolver->hullGeometry->inputVertexPositions[v] -
                                    forwardSolver->get_G()) +
                                VisualUtils::center);
    } else {
      hidden_stable_vertices.push_back(normalize(
                                           forwardSolver->hullGeometry->inputVertexPositions[v] -
                                           forwardSolver->get_G()) +
                                       VisualUtils::center);
    }
  }
  auto* stable_vertices_gm_pc =
      polyscope::registerPointCloud("stable Vertices Normals", stable_vertices);
  stable_vertices_gm_pc->setPointRadius(VisualUtils::gm_pt_radi * 3, false);
  stable_vertices_gm_pc->setPointColor({0.1, 0.9, 0.1});
  stable_vertices_gm_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

  if (VisualUtils::show_hidden_stable_vertex_normals) {
    auto* hidden =
        polyscope::registerPointCloud("hidden stable Vertices Normals", hidden_stable_vertices);
    hidden->setPointRadius(VisualUtils::gm_pt_radi, false);
    hidden->setPointColor({0.03, 0.4, 0.03});
    hidden->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  } else {
    if (polyscope::hasPointCloud("hidden stable Vertices Normals"))
      polyscope::removePointCloud("hidden stable Vertices Normals");
  }
}

void draw_stable_face_normals_on_gauss_map(Forward3DSolver* forwardSolver) {
  forwardSolver->initialize_pre_computes();
  auto* bnd_builder = new BoundaryBuilder(forwardSolver);
  bnd_builder->build_boundary_normals();

  std::vector<Vector3> stable_face_normals, face_normal_points;
  std::vector<double> stable_face_probs;
  for (Face f : forwardSolver->hullMesh->faces()) {
    Vector3 normal_pos_on_gm = forwardSolver->hullGeometry->faceNormal(f) + VisualUtils::center;
    face_normal_points.push_back(normal_pos_on_gm);
    if (forwardSolver->face_is_stable(f)) {
      stable_face_normals.push_back(normal_pos_on_gm);
      stable_face_probs.push_back(bnd_builder->face_region_area[f] / (4. * PI));
    }
  }
  auto* stable_face_normals_pc =
      polyscope::registerPointCloud("stable Face Normals", stable_face_normals);
  stable_face_normals_pc->addScalarQuantity("stable face probs", stable_face_probs);
  stable_face_normals_pc->setPointRadius(VisualUtils::gm_pt_radi * 6., false);
  stable_face_normals_pc->setPointColor({0.9, 0.1, 0.1});
  stable_face_normals_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

  auto* face_normals_pc =
      polyscope::registerPointCloud("Face Normals", face_normal_points);
  face_normals_pc->setPointRadius(VisualUtils::gm_pt_radi, false);
  face_normals_pc->setPointColor({0.9, 0.9, 0.9});
  face_normals_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  face_normals_pc->setEnabled(false);
}

void plot_height_function(Forward3DSolver* forwardSolver, ManifoldSurfaceMesh* sphere_mesh,
                          VertexPositionGeometry* sphere_geometry, bool plot_surface) {
  VertexData<Vector3> shifted_poses(*sphere_mesh);
  for (Vertex v : sphere_mesh->vertices()) {
    Vector3 pos = sphere_geometry->inputVertexPositions[v];
    shifted_poses[v] = pos.normalize() + VisualUtils::center;
  }
  VisualUtils::gm_sphere_mesh = polyscope::registerSurfaceMesh(
      "gm_sphere_mesh", shifted_poses, sphere_mesh->getFaceVertexList());
  VisualUtils::gm_sphere_mesh->setSmoothShade(true);
  VisualUtils::gm_sphere_mesh->setSurfaceColor({0.74, 0.7, 0.9});

  VertexData<double> heigh_func(*sphere_mesh);
  VertexData<Vector3> lifted_poses(*sphere_mesh);
  for (Vertex v : sphere_mesh->vertices()) {
    Vector3 pos = sphere_geometry->inputVertexPositions[v];
    double height = forwardSolver->height_function(pos.normalize());
    heigh_func[v] = height;
    lifted_poses[v] = sqrt(height) * pos.normalize() + VisualUtils::center;
  }
  auto* height_func_vis =
      VisualUtils::gm_sphere_mesh->addVertexScalarQuantity("height function", heigh_func);
  height_func_vis->setEnabled(true);

  if (plot_surface) {
    auto* height_function_mesh =
        polyscope::registerSurfaceMesh("height surface_func", lifted_poses,
                                       sphere_mesh->getFaceVertexList());
    height_function_mesh->setSmoothShade(true);
    height_function_mesh->addVertexScalarQuantity(" height function surface", heigh_func);
  }
}

void draw_gauss_map(Forward3DSolver* forwardSolver, ManifoldSurfaceMesh* sphere_mesh,
                    VertexPositionGeometry* sphere_geometry) {
  plot_height_function(forwardSolver, sphere_mesh, sphere_geometry, false);
  draw_stable_face_normals_on_gauss_map(forwardSolver);
  draw_stable_vertices_on_gauss_map(forwardSolver);
  draw_edge_arcs_on_gauss_map(forwardSolver);
  show_edge_equilibria_on_gauss_map(forwardSolver);
}

void show_edge_equilibria_on_gauss_map(Forward3DSolver* forwardSolver) {
  std::vector<Vector3> edge_equilibria_points, stabilizable_edge_equilibria_points,
      stable_edge_equilibria_points;
  Vector3 G = forwardSolver->get_G();
  for (Edge e : forwardSolver->hullMesh->edges()) {
    auto collect = [&](std::vector<Vector3>& dst) {
      Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
      Vector3 A = forwardSolver->hullGeometry->inputVertexPositions[v1],
              B = forwardSolver->hullGeometry->inputVertexPositions[v2];
      Vector3 GB = B - G, AB = B - A;
      Vector3 ortho_g = GB - AB * dot(AB, GB) / dot(AB, AB);
      dst.push_back(normalize(ortho_g) + VisualUtils::center);
    };
    if (forwardSolver->edge_is_stablizable(e) && forwardSolver->edge_is_stable(e)) {
      collect(edge_equilibria_points);
    } else if (forwardSolver->edge_is_stablizable(e)) {
      collect(stabilizable_edge_equilibria_points);
    } else if (forwardSolver->edge_is_stable(e)) {
      collect(stable_edge_equilibria_points);
    }
  }
  auto* edge_equilibria_pc =
      polyscope::registerPointCloud("Edge equilibria", edge_equilibria_points);
  edge_equilibria_pc->setPointRadius(VisualUtils::gm_pt_radi * 3., false);
  edge_equilibria_pc->setPointColor({0.2, 0.2, 0.9});
  edge_equilibria_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

  if (VisualUtils::draw_stable_g_vec_for_unstable_edge_arcs) {
    auto* stabilizable_edge_pc = polyscope::registerPointCloud(
        "Stabilizable Edge equilibria", stabilizable_edge_equilibria_points);
    stabilizable_edge_pc->setPointRadius(VisualUtils::gm_pt_radi * 3., false);
    stabilizable_edge_pc->setPointColor(VisualUtils::stabilizable_edge_color);
    stabilizable_edge_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

    auto* stable_edge_pc = polyscope::registerPointCloud("stable Edge equilibria",
                                                         stable_edge_equilibria_points);
    stable_edge_pc->setPointRadius(VisualUtils::gm_pt_radi * 3., false);
    stable_edge_pc->setPointColor(VisualUtils::stable_edge_color);
    stable_edge_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  } else {
    if (polyscope::hasPointCloud("stable Edge equilibria"))
      polyscope::removePointCloud("stable Edge equilibria");
    if (polyscope::hasPointCloud("Stabilizable Edge equilibria"))
      polyscope::removePointCloud("Stabilizable Edge equilibria");
  }
}

void draw_guess_pc(Forward3DSolver* forwardSolver,
                   std::vector<std::pair<size_t, size_t>> neigh_inds,
                   std::vector<Vector3> boundary_normals) {
  std::vector<Vector3> positions;
  for (Vector3 bnd_normal : boundary_normals) {
    double tmp_t = -1;
    for (Face f : forwardSolver->hullMesh->faces()) {
      std::vector<Vector3> face_poses;
      for (Vertex v : f.adjacentVertices())
        face_poses.push_back(forwardSolver->hullGeometry->inputVertexPositions[v]);
      tmp_t = ray_intersect(forwardSolver->get_G(), bnd_normal, face_poses);
      if (tmp_t != -1) break;
    }
    positions.push_back(forwardSolver->get_G() + tmp_t * bnd_normal);
  }
  auto* test_pc = polyscope::registerPointCloud("test point cloud", positions);
  test_pc->setPointColor({1., 0., 0.});
  test_pc->setPointRadius(VisualUtils::gm_pt_radi * 0.7, false);
  test_pc->setEnabled(true);

  std::vector<std::array<size_t, 2>> edge_inds;
  for (auto& p1 : neigh_inds) edge_inds.push_back({p1.first, p1.second});
  auto stable_graph_PSnet =
      polyscope::registerCurveNetwork("stable regions on polyhedra", positions, edge_inds);
  stable_graph_PSnet->setRadius(VisualUtils::stable_edge_radi / 2., true);
  stable_graph_PSnet->setColor(VisualUtils::stable_edge_color);
  stable_graph_PSnet->setEnabled(true);
}

void draw_G(Vector3 G) {
  VisualUtils::psG = polyscope::registerPointCloud("Center of Mass", std::vector<Vector3>{G});
  VisualUtils::psG->setPointColor({0., 0., 0.});
  VisualUtils::psG->setPointRadius(VisualUtils::G_radi, false);
}

void visualize_gauss_map(Forward3DSolver* forwardSolver) {
  std::vector<Vector3> normals_icosahedral = generate_normals_icosahedral(ICOS_RES);
  ManifoldSurfaceMesh* sphere_mesh;
  VertexPositionGeometry* sphere_geometry;
  std::tie(sphere_mesh, sphere_geometry) = get_convex_hull_mesh(normals_icosahedral);
  draw_gauss_map(forwardSolver, sphere_mesh, sphere_geometry);
}

void visualize_stable_vertices(Forward3DSolver* forwardSolver) {
  std::vector<Vector3> positions;
  for (Vertex v : forwardSolver->hullMesh->vertices()) {
    if (forwardSolver->vertex_is_stablizable(v)) {
      positions.push_back(forwardSolver->hullGeometry->inputVertexPositions[v]);
    }
  }
  auto* psCloud = polyscope::registerPointCloud("stabilizable vertices", positions);
  psCloud->setPointColor({0.1, .9, .1});
  psCloud->setPointRadius(VisualUtils::gm_pt_radi * 3., false);
  psCloud->setPointRenderMode(polyscope::PointRenderMode::Sphere);
}

void visualize_edge_stability(Forward3DSolver* forwardSolver) {
  std::vector<std::array<size_t, 2>> stable_edgeInds, stablilizable_edgeInds, both_edgeInds;
  std::vector<Vector3> stable_positions, stablilizable_positions, both_positions;
  size_t stable_counter = 0, stablizable_counter = 0, both_counter = 0;
  for (Edge e : forwardSolver->hullMesh->edges()) {
    Vector3 p1 = forwardSolver->hullGeometry->inputVertexPositions[e.firstVertex()],
            p2 = forwardSolver->hullGeometry->inputVertexPositions[e.secondVertex()];
    size_t flag = 0;
    if (forwardSolver->edge_is_stable(e)) {
      stable_positions.push_back(p1);
      stable_positions.push_back(p2);
      stable_edgeInds.push_back({stable_counter, stable_counter + 1});
      stable_counter += 2;
      flag++;
    }
    if (forwardSolver->edge_is_stablizable(e)) {
      stablilizable_positions.push_back(p1);
      stablilizable_positions.push_back(p2);
      stablilizable_edgeInds.push_back({stablizable_counter, stablizable_counter + 1});
      stablizable_counter += 2;
      flag++;
    }
    if (flag == 2) {
      both_positions.push_back(p1);
      both_positions.push_back(p2);
      both_edgeInds.push_back({both_counter, both_counter + 1});
      both_counter += 2;
    }
  }
  auto psStableEdgesNet =
      polyscope::registerCurveNetwork("stable Edges", stable_positions, stable_edgeInds);
  psStableEdgesNet->setRadius(VisualUtils::stable_edge_radi, true);
  psStableEdgesNet->setColor(VisualUtils::stable_edge_color);
  psStableEdgesNet->setEnabled(true);

  auto psStablilizableEdgesNet = polyscope::registerCurveNetwork(
      "stablizable Edges", stablilizable_positions, stablilizable_edgeInds);
  psStablilizableEdgesNet->setRadius(VisualUtils::stablizable_edge_radi, true);
  psStablilizableEdgesNet->setColor(VisualUtils::stabilizable_edge_color);
  psStablilizableEdgesNet->setEnabled(true);

  auto psBothEdgesNet =
      polyscope::registerCurveNetwork("stable && stablizable Edges", both_positions, both_edgeInds);
  psBothEdgesNet->setRadius(VisualUtils::both_edge_radi, true);
  psBothEdgesNet->setColor(VisualUtils::both_edge_color);
  psBothEdgesNet->setEnabled(true);
}

void visualize_face_stability(Forward3DSolver* forwardSolver) {
  std::vector<std::array<double, 3>> fColor(forwardSolver->hullMesh->nFaces());
  for (Face f : forwardSolver->hullMesh->faces()) {
    if (forwardSolver->face_is_stable(f))
      fColor[f.getIndex()] = {1., 0.1, 0.1};
    else
      fColor[f.getIndex()] = {0.9, 0.9, 1.};
  }
  auto* psInputMesh = polyscope::getSurfaceMesh("hull mesh");
  auto* faceQnty = psInputMesh->addFaceColorQuantity("face stability", fColor);
  faceQnty->setEnabled(true);
}

void visualize_colored_polyhedra(Forward3DSolver* forwardSolver,
                                 FaceData<Vector3> face_colors) {
  VertexData<Vector3> shifted_positions(*forwardSolver->hullMesh);
  for (Vertex v : forwardSolver->hullMesh->vertices()) {
    shifted_positions[v] =
        forwardSolver->hullGeometry->inputVertexPositions[v] + VisualUtils::colored_shift;
  }
  auto* coloredPsMesh = polyscope::registerSurfaceMesh(
      "colored polyhedra", shifted_positions, forwardSolver->hullMesh->getFaceVertexList());
  auto* faceQnty =
      coloredPsMesh->addFaceColorQuantity("random face colors", face_colors);
  faceQnty->setEnabled(true);
}

void visualize_all_stable_orientations(Forward3DSolver* forwardSolver) {
  forwardSolver->set_uniform_G();
  forwardSolver->initialize_pre_computes();
  auto* bnd_builder = new BoundaryBuilder(forwardSolver);
  bnd_builder->build_boundary_normals();

  std::vector<std::pair<Face, double>> probs;
  for (Face f : forwardSolver->hullMesh->faces())
    if (bnd_builder->face_region_area[f] > 0)
      probs.push_back({f, bnd_builder->face_region_area[f] / (4. * PI)});
  std::sort(probs.begin(), probs.end(),
            [](auto a, auto b) { return a.second > b.second; });

  double floor_z = -1.;
  Vector3 floor_vec = Vector3({0, 0, floor_z});
  VertexData<Vector3> rotated_poses(*forwardSolver->inputMesh);
  size_t i = 0;
  for (auto& p : probs) {
    Face f = p.first;
    double prob = p.second;
    Vector3 f_normal = forwardSolver->hullGeometry->faceNormal(f);
    Vector3 rot_axis = cross(f_normal, floor_vec);
    double rot_angle = angle(f_normal, floor_vec);
    Vector3 offset =
        floor_vec -
        forwardSolver->hullGeometry->inputVertexPositions[f.halfedge().vertex()]
            .rotateAround(rot_axis, rot_angle);
    for (Vertex v : forwardSolver->inputMesh->vertices()) {
      rotated_poses[v] = forwardSolver->inputGeometry->inputVertexPositions[v].rotateAround(
                             rot_axis, rot_angle) +
                         offset;
    }
    double y_shift = -(double)i;
    Vector3 vis_shift({y_shift * 2., -2., 0.});
    auto tmp_ori_mesh = polyscope::registerSurfaceMesh(
        "ori " + std::to_string(i) + " probe: " + std::to_string(prob),
        rotated_poses + vis_shift, forwardSolver->inputMesh->getFaceVertexList());
    tmp_ori_mesh->setSurfaceColor({0.1, 0.4, 0.02});
    i++;
  }
}

void draw_stable_patches_on_gauss_map(bool on_height_surface, BoundaryBuilder* bnd_builder,
                                      bool on_ambient_mesh) {
  auto net_pair =
      build_and_draw_stable_patches_on_gauss_map(bnd_builder, VisualUtils::center,
                                                 VisualUtils::gm_radi, VisualUtils::arcs_seg_count,
                                                 on_height_surface);
  if (on_ambient_mesh)
    draw_guess_pc(bnd_builder->forward_solver, net_pair.first, net_pair.second);
}

void update_visuals(Forward3DSolver* tmp_solver, BoundaryBuilder* bnd_builder,
                    ManifoldSurfaceMesh* sphere_mesh, VertexPositionGeometry* sphere_geometry) {
  draw_gauss_map(tmp_solver, sphere_mesh, sphere_geometry);
  draw_G(tmp_solver->get_G());
  plot_height_function(tmp_solver, sphere_mesh, sphere_geometry, false);
  draw_stable_patches_on_gauss_map(false, bnd_builder, false);
}

void draw_spherical_cone(std::vector<std::pair<size_t, size_t>> edges,
                         std::vector<Vector3> poses, Vector3 center, size_t seg_count,
                         glm::vec3 color, std::string title) {
  std::vector<Vector3> refined_cone_poses;
  std::vector<std::vector<size_t>> refined_cone_faces;
  refined_cone_poses.push_back(Vector3{0, 0, 0});
  size_t cone_pose_count = 1;
  for (auto& edge_pair : edges) {
    size_t i = edge_pair.first, j = edge_pair.second;
    Vector3 p1 = poses[i], p2 = poses[j];
    for (int jj = 0; jj < (int)seg_count; jj++) {
      Vector3 tmp_p1 =
                  p2 * ((double)jj / (double)seg_count) +
                  p1 * ((double)(seg_count - jj) / (double)seg_count),
              tmp_p2 =
                  p2 * ((double)(jj + 1) / (double)seg_count) +
                  p1 * ((double)(seg_count - jj - 1) / (double)seg_count);
      refined_cone_poses.push_back(tmp_p1.normalize());
      refined_cone_poses.push_back(tmp_p2.normalize());
      refined_cone_faces.push_back({0, cone_pose_count, cone_pose_count + 1});
      cone_pose_count += 2;
    }
  }
  std::vector<Vector3> refined_cone_poses_shifted;
  for (Vector3 p : refined_cone_poses) refined_cone_poses_shifted.push_back(p + center);
  polyscope::registerSurfaceMesh(title, refined_cone_poses_shifted, refined_cone_faces)
      ->setSurfaceColor(color)
      ->setBackFacePolicy(polyscope::BackFacePolicy::Identical);
}

std::tuple<ManifoldSurfaceMesh*, VertexPositionGeometry*>
make_cone_conforming_spherical_triangulation(
    BoundaryBuilder* boundary_builder, ManifoldSurfaceMesh* sphere_mesh,
    VertexPositionGeometry* sphere_geometry, Face cone_f, Vector3 shift,
    std::vector<Vector3>& cone_poses, std::vector<std::pair<size_t, size_t>>& cone_edges,
    glm::vec3 cone_color, std::vector<Vector3>& prism_poses,
    std::vector<std::pair<size_t, size_t>>& prism_edges, glm::vec3 prism_color) {
  std::vector<Vector3> s2_poses;
  for (Vertex v : sphere_mesh->vertices()) {
    Vector3 pos = sphere_geometry->inputVertexPositions[v];
    s2_poses.push_back(pos.normalize());
  }
  for (auto& edge_arc : cone_edges) {
    Vector3 p_c1 = cone_poses[edge_arc.first], p_c2 = cone_poses[edge_arc.second];
    for (Edge e : sphere_mesh->edges()) {
      Vector3 p1 = sphere_geometry->inputVertexPositions[e.firstVertex()].normalize(),
              p2 = sphere_geometry->inputVertexPositions[e.secondVertex()].normalize();
      bool sign_change = false;
      Vector3 p_int = intersect_arc_ray_with_arc(p1, p2, p_c1, p_c2, sign_change);
      if (p_int.norm() != 0 &&
          intersect_arc_ray_with_arc(p_c1, p_c2, p1, p2, sign_change).norm() != 0) {
        s2_poses.push_back(p_int);
      }
    }
  }
  for (auto& edge_arc : prism_edges) {
    Vector3 p_c1 = prism_poses[edge_arc.first], p_c2 = prism_poses[edge_arc.second];
    for (Edge e : sphere_mesh->edges()) {
      Vector3 p1 = sphere_geometry->inputVertexPositions[e.firstVertex()].normalize(),
              p2 = sphere_geometry->inputVertexPositions[e.secondVertex()].normalize();
      bool sign_change = false;
      Vector3 p_int = intersect_arc_ray_with_arc(p1, p2, p_c1, p_c2, sign_change);
      if (p_int.norm() != 0 &&
          intersect_arc_ray_with_arc(p_c1, p_c2, p1, p2, sign_change).norm() != 0) {
        s2_poses.push_back(p_int);
      }
    }
  }
  ManifoldSurfaceMesh* new_s2_mesh;
  VertexPositionGeometry* new_s2_geometry;
  std::tie(new_s2_mesh, new_s2_geometry) = get_convex_hull_mesh(s2_poses);
  polyscope::SurfaceMesh* new_s2_ps_mesh = polyscope::registerSurfaceMesh(
      "new_s2_mesh", new_s2_geometry->inputVertexPositions + shift,
      new_s2_mesh->getFaceVertexList());
  polyscope::getSurfaceMesh("gm_sphere_mesh")->setEnabled(false);
  new_s2_ps_mesh->setSmoothShade(true);
  new_s2_ps_mesh->setSurfaceColor({0.74, 0.7, 0.9});
  new_s2_ps_mesh->setTransparency(0.5);

  FaceData<Vector3> face_colors(*new_s2_mesh, Vector3({0.74, 0.7, 0.9}));
  Vector3 face_cone_p1 = prism_poses[1].normalize(), face_cone_p2 = prism_poses[2].normalize(),
          face_cone_p3 = prism_poses[3].normalize();
  std::vector<std::vector<size_t>> face_cone_faces, MS_cone_faces;
  std::vector<std::vector<size_t>> all_faces_list = new_s2_mesh->getFaceVertexList();
  for (Face f : new_s2_mesh->faces()) {
    Vector3 p1 = new_s2_geometry->inputVertexPositions[f.halfedge().vertex()].normalize(),
            p2 = new_s2_geometry->inputVertexPositions[f.halfedge().next().vertex()].normalize(),
            p3 = new_s2_geometry->inputVertexPositions[f.halfedge().next().next().vertex()]
                     .normalize();
    Vector3 center_normal = ((p1 + p2 + p3) / 3.).normalize();
    Face final_face = boundary_builder->forward_solver->final_touching_face(center_normal);
    if (final_face == cone_f) {
      face_colors[f] = Vector3({cone_color.x, cone_color.y, cone_color.z});
      MS_cone_faces.push_back(all_faces_list[f.getIndex()]);
    }
    if (is_in_positive_cone(face_cone_p1, face_cone_p2, face_cone_p3, center_normal)) {
      face_colors[f] = Vector3({prism_color.x, prism_color.y, prism_color.z});
      face_cone_faces.push_back(all_faces_list[f.getIndex()]);
    }
  }

  new_s2_ps_mesh->addFaceColorQuantity("cone MS colors", face_colors)->setEnabled(true);
  polyscope::registerSurfaceMesh("MS cone patch", new_s2_geometry->inputVertexPositions + shift,
                                 MS_cone_faces)
      ->setSurfaceColor(cone_color);
  polyscope::registerSurfaceMesh("Face solid angle patch",
                                 new_s2_geometry->inputVertexPositions + shift, face_cone_faces)
      ->setSurfaceColor(prism_color);
  return {new_s2_mesh, new_s2_geometry};
}

void visualize_face_solid_angle_vs_ms_complex(size_t f_ind, BoundaryBuilder* boundary_builder,
                                              ManifoldSurfaceMesh* sphere_mesh,
                                              VertexPositionGeometry* sphere_geometry) {
  Face f = boundary_builder->forward_solver->hullMesh->face(f_ind);
  if (boundary_builder->forward_solver->face_last_face[f] != f)
    throw std::logic_error("provide stable face for plot\n");

  Vector3 shift({0, 0, 0});
  std::vector<Vector3> prism_verts;
  Vector3 G = boundary_builder->forward_solver->get_G();
  prism_verts.push_back(G);
  for (Vertex v : f.adjacentVertices()) {
    Vector3 pos = boundary_builder->forward_solver->hullGeometry->inputVertexPositions[v];
    prism_verts.push_back(pos.normalize());
  }
  std::vector<std::pair<size_t, size_t>> prism_edge_pairs;
  prism_edge_pairs.push_back({1, 2});
  prism_edge_pairs.push_back({2, 3});
  prism_edge_pairs.push_back({3, 1});
  glm::vec3 prism_color = glm::vec3({1., 14. / 25., 0.});
  draw_arc_network_on_sphere(prism_edge_pairs, prism_verts, shift, 1., 12,
                             "Face Solid angle cone", 1., prism_color, 0.009);
  draw_spherical_cone(prism_edge_pairs, prism_verts, shift, 12, prism_color,
                      "refined face solid angle cone");

  auto p = boundary_builder->MS_complex_edges_of_face(f);
  std::vector<Vector3> boundary_normals = p.second;
  std::vector<std::pair<size_t, size_t>> ind_pairs = p.first;
  glm::vec3 MS_cone_color = glm::vec3({33. / 255., 100. / 255., 10. / 255.});
  draw_arc_network_on_sphere(ind_pairs, boundary_normals, shift, 1., 12, "MS Cone boundary", 1.,
                             MS_cone_color, 0.009);
  draw_spherical_cone(ind_pairs, boundary_normals, shift, 12, MS_cone_color, "refined MS cone");

  make_cone_conforming_spherical_triangulation(boundary_builder, sphere_mesh, sphere_geometry, f,
                                               shift, boundary_normals, ind_pairs, MS_cone_color,
                                               prism_verts, prism_edge_pairs,
                                               prism_color);
}


void init_visuals(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry,
				  Forward3DSolver* forwardSolver, BoundaryBuilder* boundary_builder){
	auto psInputMesh = polyscope::registerSurfaceMesh(
		"init input mesh",
		geometry->inputVertexPositions, mesh->getFaceVertexList())->setTransparency(1.0);
	auto psHullMesh = polyscope::registerSurfaceMesh(
		"init hull mesh", forwardSolver->hullGeometry->inputVertexPositions, 
		forwardSolver->hullMesh->getFaceVertexList())->setEdgeWidth(0.7)->setTransparency(0.9);
	draw_G(forwardSolver->get_G());
	visualize_gauss_map(forwardSolver);
  draw_stable_patches_on_gauss_map(boundary_builder, false);

}


void draw_stable_patches_on_gauss_map(BoundaryBuilder* boundary_builder,
									  bool on_height_surface){
  auto net_pair = build_and_draw_stable_patches_on_gauss_map(boundary_builder,
                                                              VisualUtils::center, VisualUtils::gm_radi, VisualUtils::arcs_seg_count, 
                                                              on_height_surface);
}


void draw_trail_on_gm(std::vector<Vector3> trail, glm::vec3 color, std::string name, double radi, bool color_gradient){
  std::vector<std::pair<size_t, size_t>> edge_inds;
  if (!color_gradient){
    for (size_t i = 0; i < trail.size()-1; i++)
      edge_inds.push_back({i, i+1});
    draw_arc_network_on_sphere(edge_inds, trail, VisualUtils::center, VisualUtils::gm_radi, VisualUtils::arcs_seg_count, name, radi, color);//{0.7,0.1,0.8}
  }
  else{
    std::vector<glm::vec3> colors;
    for (size_t i = 0; i < trail.size(); i++){
      glm::vec3 tmp_color = color;
      tmp_color.x = tmp_color.x * ((i/(float)trail.size()));
      draw_arc_on_sphere(trail[i], trail[i+1], VisualUtils::center, VisualUtils::gm_radi, VisualUtils::arcs_seg_count, i, radi, tmp_color);
    }
  }
}


void visualize_quasi_static_drop_sequence(
                    std::vector<Eigen::Matrix4d> transformation_matrices,
					std::vector<Vector3> saved_snail_trail_refined,
					Forward3DSolver* forwardSolver
					){
	std::cout << "visualizing snail trail steps \n";
	VertexData<Vector3> tmp_positions(*forwardSolver->inputMesh);
	// tmp_positions = forwardSolver->inputGeometry->inputVertexPositions;
	for (int i = 0; i < transformation_matrices.size(); i++){
		Eigen::Matrix4d transformation_matrix = transformation_matrices[i];
		// apply to the input mesh
		for (Vertex v: forwardSolver->inputMesh->vertices()){
			Vector3 p = forwardSolver->inputGeometry->inputVertexPositions[v];
			// apply the transformation matrix
			Eigen::Vector4d tmp_v({p.x, p.y, p.z, 1.0});
			tmp_v = transformation_matrix * tmp_v;
			tmp_positions[v] = Vector3({tmp_v[0], tmp_v[1], tmp_v[2]});
		}
		// visualize
		polyscope::registerSurfaceMesh("quasi static snail step ", 
			tmp_positions, 
			forwardSolver->inputMesh->getFaceVertexList())->setEnabled(true);
		// only the snail trail
		Vector3 tmp_orientation = saved_snail_trail_refined[i];
		std::vector<Vector3> tmp_orientation_vec = {tmp_orientation * VisualUtils::gm_radi + VisualUtils::center};
		polyscope::registerPointCloud("quasi orientation on gm", 
										tmp_orientation_vec)->setPointRadius(0.03, false)->setPointColor({0,0,0});
		polyscope::show();
	}  
}

void draw_ground_plane_mesh(Vector3 down_vec, double height, double half_width){
  std::vector<std::vector<size_t>> faces = {{0, 1, 2, 3}}; // single quad
  Vector3 ux{1, 0, 0}, uy;
  if (cross(down_vec, ux).norm() < 1e-6) {
    ux = Vector3({0, 1, 0}); 
  }
  ux = ux - dot(ux, down_vec) * down_vec; // project ux onto the plane
  ux = ux.normalize();
  uy = cross(down_vec, ux).normalize(); // orthogonal to both down_vec and ux

  std::vector<Vector3> positions = {
      down_vec * height + half_width * ux + half_width * uy, // 1,1
      down_vec * height + half_width * ux - half_width * uy, // 1,-1
      down_vec * height - half_width * ux - half_width * uy, // -1,-1
      down_vec * height - half_width * ux + half_width * uy  // -1,1
  };
  polyscope::registerSurfaceMesh("ground plane", positions, faces)->setSurfaceColor({0.5, 0.8, 0.5})->setTransparency(0.7);
}