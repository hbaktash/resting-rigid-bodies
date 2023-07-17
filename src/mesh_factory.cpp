#include "mesh_factory.h"

Vector3 spherical_to_xyz(double r, double phi, double theta){
  return Vector3({r*cos(phi)*sin(theta), r*cos(phi)*cos(theta), r*sin(phi)});
}





std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> 
generate_polyhedra(std::string poly_str){
    std::vector<std::vector<size_t>> faces;
    std::vector<Vector3> positions;
    int n;
    if (std::strcmp(poly_str.c_str(), "tet") == 0){
      n = 4;
      faces = {{0, 2, 1},
                              {0, 1, 3},
                              {0, 3, 2},
                              {2, 3, 1}};
      double theta0 = 0., theta1 = 2.*PI/3., theta2 = 4.*PI/3.,
             phi0 = PI/6.,
             theta3 = 0., 
             phi1 = -PI/2.;
      
      positions = { spherical_to_xyz(1., phi0, theta0),
                    spherical_to_xyz(1., phi0, theta1),
                    spherical_to_xyz(1., phi0, theta2),
                    spherical_to_xyz(1., phi1, theta3)};
    }
    else if (std::strcmp(poly_str.c_str(), "tet2") == 0){
      n = 4;
      faces = {{0, 1, 2},
                              {0, 3, 1},
                              {0, 2, 3},
                              {2, 1, 3}};
      double theta0 = 0., theta1 = 2.*PI/3., theta2 = 4.*PI/3.,
             phi0 = PI/6.,
             theta3 = 0., 
             phi1 = -PI/2.;
      
      positions = { Vector3({0.1 , -1.5, -0.2})/3.,
                    Vector3({-1., 1., 1.2})/3.,
                    Vector3({-1.1, 1., -1.})/3.,
                    Vector3({2. , 1.5, 0.})/3.};
    }
    else if (std::strcmp(poly_str.c_str(), "sliced tet") == 0){
      double theta0 = 0., theta1 = 2.*PI/3., theta2 = 4.*PI/3.,
              phi0 = PI/6.,
              theta3 = 0., 
              phi1 = -PI/2.;
      Vector3 p1 = spherical_to_xyz(1., phi0, theta0), //0
              p2 = spherical_to_xyz(1., phi0, theta1), // gone
              p3 = spherical_to_xyz(1., phi0, theta2), // gone
              p = spherical_to_xyz(1., phi1, theta3); // 1
      Vector3 p2p  = 0.8*p2 + 0.2*p, // 2
              p2p3 = 0.7*p2 + 0.3*p3, // 3
              p2p1 = 0.8*p2 + 0.2*p1; // 4
      Vector3 p3p  = 0.8*p3 + 0.2*p, // 5
              p3p2 = 0.7*p3 + 0.3*p2, // 6
              p3p1 = 0.8*p3 + 0.2*p1; // 7
      // making the quad faces non-planar
      double shrinkage = 0.7;
      p2p *= 0.8; 
      p2p1 *= shrinkage;
      p3p *= 0.8; 
      p3p1 *= shrinkage;
      positions.push_back(p); positions.push_back(p1);
      positions.push_back(p2p); positions.push_back(p2p3); positions.push_back(p2p1); 
      positions.push_back(p3p); positions.push_back(p3p2); positions.push_back(p3p1);
      faces = {{0, 3, 6},
               {0, 2, 3},
               {0, 6, 5},
               {0, 1, 2},
               {0, 5, 1},
               {1, 3, 4},
               {1, 7, 6},   
               {1, 6, 3},
               {5, 6, 7},
               {2, 4, 3},
               {1, 4, 2},
               {1, 5, 7}};
    }
    else if (std::strcmp(poly_str.c_str(), "rndbs 6-gon 1") == 0){
    }
    else if (std::strcmp(poly_str.c_str(), "rndbs 9-gon 1") == 0){
    }
    else {
      throw std::runtime_error("no valid string provided\n");
    }
    // std::unique_ptr<ManifoldSurfaceMesh> poly_triangulated;
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    mesh.reset(new ManifoldSurfaceMesh(faces));
    std::unique_ptr<VertexPositionGeometry> geometry(new VertexPositionGeometry(*mesh));
    for (Vertex v : mesh->vertices()) {
        // Use the low-level indexers here since we're constructing
        (*geometry).inputVertexPositions[v] = positions[v.getIndex()];
    }
    return std::make_tuple(std::move(mesh), std::move(geometry));
}
