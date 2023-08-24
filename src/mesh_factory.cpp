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
#include "mesh_factory.h"


Vector3 spherical_to_xyz(double r, double phi, double theta){
  return Vector3({r*cos(phi)*sin(theta), r*cos(phi)*cos(theta), r*sin(phi)});
}

Vector3 cylindrical_to_xyz(double h, double r, double theta){        
  return Vector3({r*cos(theta), h, r*sin(theta)});
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
    else if (std::strcmp(poly_str.c_str(), "Conway spiral 4") == 0){
      // from https://arxiv.org/abs/2103.13727
      std::vector<double> angles({66.173, 44.519, 29.875, 19.716, 19.716});
      size_t tip_ind = 5 * 4 + 1 - 1;// 21 vertices total
      double tmp_r = 1.,
             tmp_theta = 0.;
      for (size_t i = 0; i < 4; i++){ // the 2d spiral
        double tmp_alpha = angles[i]*PI/180.;
        tmp_theta += tmp_alpha;
        printf("tmp theta %f\n", tmp_theta);
        tmp_r *= cos(tmp_alpha);
        double cyli_r = tmp_r * sin(tmp_theta),
               cyli_h = tmp_r * cos(tmp_theta);
        for (size_t j = 0; j < 5; j++){ // rotation along an axis by PI/5
          double cyli_theta = (double)j * 2.* PI/ 5.;
          positions.push_back(cylindrical_to_xyz(cyli_h, cyli_r, cyli_theta));
          size_t ind = 5 * i + j,
                 next_ind = 5 * i + ((j+1) % 5); 
          if (i == 0){
            faces.push_back({ind, tip_ind, next_ind});
          }
          else {
            faces.push_back({next_ind, ind, ind - 5, next_ind - 5});
          }
        }
      }
      // bottom face
      faces.push_back({15,16,17,18,19});
      // faces.push_back({19,18,17,16,15});
      // tip vertex location
      positions.push_back(cylindrical_to_xyz(1., 0., 0.));

      for (std::vector<size_t> f: faces){
        for (size_t ind: f){
          printf(" %d,", ind);
        }
        printf("\n");
      }
      for (Vector3 pos: positions){
        std::cout << "pos: " << pos << "\n";
      }
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
