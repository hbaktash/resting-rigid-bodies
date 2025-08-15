#include "inverse_design/prob_assignment.h"



std::vector<std::pair<Vector3, double>> 
normal_prob_assignment(std::string shape_name){
  std::vector<std::pair<Vector3, double>> normal_to_prob_pairs;
  
  if (shape_name == "circus"){
    normal_to_prob_pairs = {
      {Vector3({0,0,-1})              , 6./36.},
      {Vector3({0.629,0.457,0.629})   , 5./36.},
      {Vector3({0.809,0.587,0.})      , 5./36.},
      {Vector3({-0.24,0.739,0.629})   , 4./36.},
      {Vector3({-0.309, 0.951, 0.})   , 4./36.},
      {Vector3({-0.777, 0., 0.629})   , 3./36.},
      {Vector3({-1., 0., 0.})         , 3./36.},
      {Vector3({-0.24, -0.739, 0.629}), 2./36.},
      {Vector3({-0.309, -0.951, 0.})  , 2./36.},
      {Vector3({0.629, -0.456, 0.629}), 1./36.},
      {Vector3({0.809, -0.587, 0.})   , 1./36.}
    };
  }
  else if (shape_name == "hendecahedron"){
    normal_to_prob_pairs = {
      // {Vector3({0, 0, -1})                    , 6./36.}, // top unique quad
      // {Vector3({0, -0.707107, 0.707107})      , 5./36.}, // bottom 4 faces
      // {Vector3({0, 0.707107, 0.707107})       , 5./36.}, 
      // {Vector3({-0.57735, -0.57735, -0.57735}), 4./36.}, // top, right quads
      // {Vector3({0.57735, -0.57735, -0.57735}) , 4./36.}, 
      // {Vector3({-0.57735, 0.57735, -0.57735}) , 3./36.}, // top, left quads
      // {Vector3({0.57735, 0.57735, -0.57735})  , 3./36.}, 
      // {Vector3({-0.894427, 0, -0.447214})     , 2./36.}, // top ring; 6 faces
      // {Vector3({0.894427, 0, -0.447214})      , 2./36.}, // triangle   
      // {Vector3({0.447214, 0, 0.894427})       , 1./36.}, // triangle
      // {Vector3({-0.447214, 0, 0.894427})      , 1./36.}  // triangle
      // 0.16: 6/36, 0.13: 5/36, 0.11: 4/36, 0.08: 3/36, 0.05: 2/36, 0.03: 1/36
      {Vector3({0, 0, -1})                    , 6./36.}, // top unique quad
      {Vector3({-0.57735, -0.57735, -0.57735}), 4./36.}, // top, right quads
      {Vector3({0.57735, -0.57735, -0.57735}) , 4./36.}, 
      {Vector3({-0.57735, 0.57735, -0.57735}) , 5./36.}, // top, left quads
      {Vector3({0.57735, 0.57735, -0.57735})  , 5./36.}, 
      {Vector3({-0.894427, 0, -0.447214})     , 3./36.}, // top ring; 6 faces
      {Vector3({0.894427, 0, -0.447214})      , 3./36.}, // triangle   
      {Vector3({0, -0.707107, 0.707107})      , 2./36.}, // bottom 4 faces
      {Vector3({0, 0.707107, 0.707107})       , 2./36.}, 
      {Vector3({0.447214, 0, 0.894427})       , 1./36.}, // triangle
      {Vector3({-0.447214, 0, 0.894427})      , 1./36.}  // triangle
    };
  }
  else if (shape_name == "cube"){
    normal_to_prob_pairs = {
      {Vector3({0,  0 ,-1}), binomial_dist(5, 0)},
      {Vector3({0,  0 , 1}), binomial_dist(5, 5)},
      {Vector3({0, -1,  0}), binomial_dist(5, 4)},
      {Vector3({0,  1,  0}), binomial_dist(5, 1)},
      {Vector3({ 1, 0,  0}), binomial_dist(5, 2)},
      {Vector3({-1, 0,  0}), binomial_dist(5, 3)}
    };
  }
  else if (shape_name == "dodecahedron"){
    normal_to_prob_pairs = {
      // antipodal formation
      {Vector3({-0.525731,  0.850651, 0}),  binomial_dist(11, 5)},
      {Vector3({ 0.525731, -0.850651, 0}),  binomial_dist(11, 5)},
      {Vector3({ 0.850651, 0, -0.525731}),  binomial_dist(11, 4)},
      {Vector3({-0.850651, 0,  0.525731}),  binomial_dist(11, 4)},
      {Vector3({ 0, -0.525731, -0.850651}), binomial_dist(11, 0)},
      {Vector3({ 0,  0.525731,  0.850651}), binomial_dist(11, 0)},
      {Vector3({ 0.525731,  0.850651, 0}),  binomial_dist(11, 1)},
      {Vector3({-0.525731, -0.850651, 0}),  binomial_dist(11, 1)},
      {Vector3({ 0,  0.525731, -0.850651}), binomial_dist(11, 3)},
      {Vector3({ 0, -0.525731,  0.850651}), binomial_dist(11, 3)},
      {Vector3({ 0.850651, 0,  0.525731}),  binomial_dist(11, 2)},
      {Vector3({-0.850651, 0, -0.525731}),  binomial_dist(11, 2)}
      // DwarfRing formation
      // {Vector3({ 0.525731, -0.850651, 0}),  binomial_dist(11, 5)},
      // {Vector3({ 0, -0.525731, -0.850651}), binomial_dist(11, 6)},
      // {Vector3({ 0.850651, 0, -0.525731}),  binomial_dist(11, 4)},
      // {Vector3({-0.525731, -0.850651, 0}),  binomial_dist(11, 7)},

      // {Vector3({ 0.850651, 0,  0.525731}),  binomial_dist(11, 3)},
      // {Vector3({ 0, -0.525731,  0.850651}), binomial_dist(11, 8)},

      // {Vector3({-0.525731,  0.850651, 0}),  binomial_dist(11, 0)},
      // {Vector3({ 0,  0.525731,  0.850651}), binomial_dist(11, 11)},

      // {Vector3({ 0.525731,  0.850651, 0}),  binomial_dist(11, 1)},
      // {Vector3({-0.850651, 0,  0.525731}),  binomial_dist(11, 10)},

      // {Vector3({-0.850651, 0, -0.525731}),  binomial_dist(11, 2)},
      // {Vector3({ 0,  0.525731, -0.850651}), binomial_dist(11, 9)},
    };
  }
  else if (shape_name == "octahedron"){
    normal_to_prob_pairs = {
      {Vector3({ 0.816497, 0, -0.57735}),  binomial_dist(7, 3)},
      {Vector3({-0.816497, 0, -0.57735}),  binomial_dist(7, 2)},
      {Vector3({0, -0.816497, -0.57735}),  binomial_dist(7, 4)},
      {Vector3({0,  0.816497, -0.57735}),  binomial_dist(7, 5)},
      {Vector3({0, -0.816497,  0.57735}),  binomial_dist(7, 1)},
      {Vector3({0,  0.816497,  0.57735}),  binomial_dist(7, 0)},
      {Vector3({ 0.816497, 0,  0.57735}),  binomial_dist(7, 6)},
      {Vector3({-0.816497, 0,  0.57735}),  binomial_dist(7, 7)},
    };
  }
  else if (shape_name.substr(shape_name.find(" ") + 1) == "prism"){
    size_t sides = std::stoi(shape_name.substr(0, shape_name.find(" ")));
    for (size_t i = 0; i < sides; i++){
      double angle1 = 2.*PI*(double)i/(double)sides,
             angle2 = 2.*PI*(double)((i + 1) % sides)/(double)sides;
      Vector3 p1 = Vector3({cos(angle1), sin(angle1), 0}),
              p2 = Vector3({cos(angle2), sin(angle2), 0});
      Vector3 normal = (p1 + p2).normalize();
      
      int binom_k = -1,
          n = sides - 1;
      // antipodal assignment
      // if (i <= n/2)
      //   binom_k = i;
      // else
      //   binom_k = i - n/2 - 1;
      binom_k = i; //  wide approach
      normal_to_prob_pairs.push_back({normal, binomial_dist(sides - 1, binom_k)});
    }
  }

  return normal_to_prob_pairs;

  // // wide tent
  // Vector3  nf7({0, 0, -1}); // bottom face
  // Vector3 nf61({0.803215, 0.26098, 0.535477});   // 0
  // Vector3 nf62({0.803215, -0.26098, 0.535477});  // 9
  // Vector3 nf51({0.496414, 0.683255, 0.535477});  // 1
  // Vector3 nf52({0.496414, -0.683255, 0.535477}); // 8
  // Vector3 nf41({0, 0.84455, 0.535477});          // 2
  // Vector3 nf42({0, -0.84455, 0.535477});         // 7
  // Vector3 nf31({-0.496414, 0.683255, 0.535477}); // 3
  // Vector3 nf32({-0.496414, -0.683255, 0.535477});// 6
  // Vector3 nf21({-0.803215, 0.26098, 0.535477});  // 4
  // Vector3 nf22({-0.803215, -0.26098, 0.535477}); // 5

  // atipodal tent
  // Vector3 nf7({0, 0, -1}); // bottom face
  // Vector3 nf61({0.803215, 0.26098, 0.535477});   // 0
  // Vector3 nf62({-0.803215, -0.26098, 0.535477}); // 5
  // Vector3 nf51({0.496414, 0.683255, 0.535477});  // 1
  // Vector3 nf52({-0.496414, -0.683255, 0.535477});// 6
  // Vector3 nf31({0, 0.84455, 0.535477});          // 2
  // Vector3 nf32({0, -0.84455, 0.535477});         // 7
  // Vector3 nf21({-0.496414, 0.683255, 0.535477}); // 3
  // Vector3 nf22({0.496414, -0.683255, 0.535477}); // 8
  // Vector3 nf41({-0.803215, 0.26098, 0.535477});  // 4
  // Vector3 nf42({0.803215, -0.26098, 0.535477});  // 9
}


// initializer for fair cluster assignment 
std::vector<std::pair<Vector3, double>> 
normal_prob_assignment_fair(Forward3DSolver *tmp_solver, size_t dice_side_count){
  BoundaryBuilder tmp_bnd_builder(tmp_solver);
  tmp_bnd_builder.build_boundary_normals();
  // sort face region areas
  std::vector<std::pair<Face, double>> face_areas;
  for (Face f: tmp_solver->hullMesh->faces()){
    if (tmp_bnd_builder.face_region_area[f] > 0){
      face_areas.push_back({f, tmp_bnd_builder.face_region_area[f]});
    }
  }
  std::sort(face_areas.begin(), face_areas.end(), [] (auto a, auto b) { return a.second > b.second; });

  // assign
  std::vector<std::pair<Vector3, double>> normal_prob_pairs;
  double goal_prob = 1./(double)dice_side_count;
  for (size_t i = 0; i < dice_side_count; i++){
    if (i >= face_areas.size()){
      break;
    }
    Face f = face_areas[i].first;
    Vector3 normal = tmp_solver->hullGeometry->faceNormal(f);
    normal_prob_pairs.push_back({normal, goal_prob});
  }
  return normal_prob_pairs;
}


FaceData<double> 
manual_stable_only_face_prob_assignment(Forward3DSolver *tmp_solver, std::vector<std::pair<Vector3, double>> normal_prob_pairs){
  FaceData<double> goal_probs(*tmp_solver->hullMesh, 0.);

  std::vector<std::pair<Vector3, double>> normal_to_prob_pairs = normal_prob_pairs;
  
  for (auto normal_prob_pair: normal_to_prob_pairs){
    Vector3 normal = normal_prob_pair.first;
    double goal_prob = normal_prob_pair.second;
    double smallest_normal_diff = 1000.;
    Face closest_face;
    for (Face f: tmp_solver->hullMesh->faces()){
      double normal_diff = (tmp_solver->hullGeometry->faceNormal(f) - normal).norm();
      if (tmp_solver->face_is_stable(f) && 
          normal_diff < smallest_normal_diff){
        closest_face = f;
        smallest_normal_diff = normal_diff;
      }
    }
    goal_probs[closest_face] = goal_prob;
  }
  return goal_probs;
}

std::vector<std::tuple<std::vector<Face>, double, Vector3>> 
manual_clustered_face_prob_assignment(Forward3DSolver *tmp_solver, std::vector<std::pair<Vector3, double>> normal_prob_pairs){
  // std::string shape_name = "octahedron binomial"; // "circus", "hendecahedron", "wide tent", "atipodal tent", "icosahedron binomial", "cube binomial", dodecahedron binomial
  
  FaceData<Vector3> closest_normals(*tmp_solver->hullMesh);
  for (Face f: tmp_solver->hullMesh->faces()){
    Vector3 closest_normal;
    double smallest_normal_diff = 1000.;
    for (auto normal_prob_pair: normal_prob_pairs){ // assign a probable normal to each face (stable/unstable)
      Vector3 normal = normal_prob_pair.first;
      double normal_diff = (tmp_solver->hullGeometry->faceNormal(f) - normal).norm();
      if (normal_diff < smallest_normal_diff){
        closest_normal = normal;
        smallest_normal_diff = normal_diff;
      }
    }
    closest_normals[f] = closest_normal;
  }

  std::vector<std::tuple<std::vector<Face>, double, Vector3>> clustered_prob_assignment; // to be populated and returned
  for (auto normal_prob_pair: normal_prob_pairs){
    
    Vector3 normal = normal_prob_pair.first;
    std::vector<Face> assigned_faces;
    double goal_prob = normal_prob_pair.second;

    for (Face f: tmp_solver->hullMesh->faces()){
      if (closest_normals[f] == normal){
        assigned_faces.push_back(f);
      }
    }
    clustered_prob_assignment.push_back(std::make_tuple(assigned_faces, goal_prob, normal));
  }
  return clustered_prob_assignment;
}

std::vector<std::pair<Vector3, double>> 
update_normal_prob_assignment(Forward3DSolver *tmp_solver,
                              std::vector<std::tuple<std::vector<Face>, double, Vector3>> clustered_face_normals,
                              bool take_max_prob_face){
  std::vector<std::pair<Vector3, double>> updated_assignment;
  for (auto cluster: clustered_face_normals){
    std::vector<Face> faces = std::get<0>(cluster);
    double cluster_goal_prob = std::get<1>(cluster);
    Vector3 normal = std::get<2>(cluster);
    
    BoundaryBuilder tmp_bnd_builder(tmp_solver);
    tmp_bnd_builder.build_boundary_normals();
    // // Two choices
    if (!take_max_prob_face){
      Vector3 cluster_stables_avg{0., 0., 0.};
      double prob_sum = 0.;
      // average the faces
      for (Face f: faces){
        if (tmp_solver->face_last_face[f] == f){
          cluster_stables_avg += tmp_solver->hullGeometry->faceNormal(f) * tmp_bnd_builder.face_region_area[f];
          prob_sum += tmp_bnd_builder.face_region_area[f];
        }
      }
      cluster_stables_avg /= prob_sum;
      cluster_stables_avg = cluster_stables_avg.normalize();
      updated_assignment.push_back({cluster_stables_avg, cluster_goal_prob});
    }
    else{
      // take max prob face
      Face max_prob_face;
      double max_prob = 0.;
      for (Face f: faces){
        if (tmp_bnd_builder.face_region_area[f] > max_prob){
          max_prob_face = f;
          max_prob = tmp_bnd_builder.face_region_area[f];
        }
      }
      updated_assignment.push_back({tmp_solver->hullGeometry->faceNormal(max_prob_face), cluster_goal_prob});
    }

  }
  return updated_assignment;
}

