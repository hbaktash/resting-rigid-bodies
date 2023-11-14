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
#include "geometrycentral/surface/meshio.h"


Vector3 spherical_to_xyz(double r, double phi, double theta){
  return Vector3({r*cos(phi)*sin(theta), r*cos(phi)*cos(theta), r*sin(phi)});
}

Vector3 cylindrical_to_xyz(double h, double r, double theta){        
  return Vector3({r*cos(theta), h, r*sin(theta)});
}

// for oloid; reindex to avoid duplicate vertices (otherwise, mesh will have a boundary)
int idx( int i, int n ) {
   if( i > 0 )
      i--;

   if( i == n*4-2 )
      i--;

   return i;
}

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> 
generate_polyhedra(std::string poly_str){
    std::vector<std::vector<size_t>> faces;
    std::vector<Vector3> positions;
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;
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
      faces = {{0, 2, 1},
                              {0, 1, 3},
                              {0, 3, 2},
                              {2, 3, 1}};
      double theta0 = 0., theta1 = 1.*PI/3., theta2 = 4.*PI/3.,
             phi0 = PI/6.,
             theta3 = 0., 
             phi1 = -PI/2. + PI/8.;
      
      positions = { spherical_to_xyz(1., phi0, theta0),
                    spherical_to_xyz(1., phi0, theta1),
                    spherical_to_xyz(1., phi0, theta2),
                    spherical_to_xyz(1., phi1, theta3)};
    }
    else if (std::strcmp(poly_str.c_str(), "cube") == 0){
      
      // bottom face
      positions.push_back(Vector3({-1, -1, -1})/2.);
      positions.push_back(Vector3({-1,  1, -1})/2.);
      positions.push_back(Vector3({ 1,  1, -1})/2.);
      positions.push_back(Vector3({ 1, -1, -1})/2.);
      // top face
      positions.push_back(Vector3({-1, -1,  1})/2.);
      positions.push_back(Vector3({-1,  1,  1})/2.);
      positions.push_back(Vector3({ 1,  1,  1})/2.);
      positions.push_back(Vector3({ 1, -1,  1})/2.);
      // v0 adjs
      faces.push_back({0, 1, 2, 3});
      faces.push_back({0, 4, 5, 1});
      faces.push_back({0, 3, 7, 4});
      // v6 adjs
      faces.push_back({6, 5, 4, 7});
      faces.push_back({6, 2, 1, 5});
      faces.push_back({6, 7, 3, 2});
    }
    else if (std::strcmp(poly_str.c_str(), "tilted cube") == 0){
      // bottom face
      positions.push_back(Vector3({-1, -1, -1})/4.);
      positions.push_back(Vector3({-1,  1, -1})/4.);
      positions.push_back(Vector3({ 1,  1, -1})/4.);
      positions.push_back(Vector3({ 1, -1, -1})/4.);
      // top face
      positions.push_back(Vector3({-1, -1,  1})/2.);
      positions.push_back(Vector3({-1,  1,  1})/2.);
      positions.push_back(Vector3({ 1,  1,  1})/2.);
      positions.push_back(Vector3({ 1, -1,  1})/2.);
      // v0 adjs
      faces.push_back({0, 1, 2, 3});
      faces.push_back({0, 4, 5, 1});
      faces.push_back({0, 3, 7, 4});
      // v6 adjs
      faces.push_back({6, 5, 4, 7});
      faces.push_back({6, 2, 1, 5});
      faces.push_back({6, 7, 3, 2});
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

      // for (std::vector<size_t> f: faces){
      //   for (size_t ind: f){
      //     printf(" %d,", ind);
      //   }
      //   printf("\n");
      // }
      // for (Vector3 pos: positions){
      //   std::cout << "pos: " << pos << "\n";
      // }
    }
    else if (std::strcmp(poly_str.c_str(), "oloid") == 0){
      int nPolygons = 25;
      int N = nPolygons/4 + 1;

      // output vertices
      for( int i = 0; i < N; i++ ) {
        const double pi = 3.141592653589793;
        double s = (double) i / (double) (N-1);
        double t = (2.*pi/3.) * ( 3.*s - 2.*pow( s, 1.5 ) ); // reparameterize to get more uniform spacing

        // circular arc 1
        double ax = sin(t);
        double ay = -cos(t);
        double az = 0.;

        // circular arc 2
        double bx = 0.;
        double by = 1./(1.+cos(t));
        double bz = sqrt(1.+2.*cos(t))/(1.+cos(t));

        // construct each arc using two separate pieces (for symmetry)
        positions.push_back(Vector3({ax,ay,az}));
        if(i > 0) positions.push_back(Vector3({-ax,ay,az}));
        positions.push_back(Vector3({bx,by,bz}));
        if(i < N - 1) positions.push_back(Vector3({bx,by,-bz}));
        //               out << "v " <<  ax << " " << ay << " " <<  az << endl;
        // if( i  >  0  ) out << "v " << -ax << " " << ay << " " <<  az << endl; // omit duplicate vertex
        //               out << "v " <<  bx << " " << by << " " <<  bz << endl;
        // if( i != N-1 ) out << "v " <<  bx << " " << by << " " << -bz << endl; // omit duplicate vertex
      }

      // output quads, as pairwise connections between each of two pieces on two arcs
      for( int i = 0; i < N-1; i++ ) {
        size_t v0, v1, v2, v3;
        size_t j = i+1;

        v0 = idx( j*4+0, N );
        v1 = idx( j*4+2, N );
        v2 = idx( i*4+2, N );
        v3 = idx( i*4+0, N );
        faces.push_back({v0, v1, v2, v3});
        // out << "f " << v0 << " " << v1 << " " << v2 << " " << v3 << endl;

        v0 = idx( i*4+0, N );
        v1 = idx( i*4+3, N );
        v2 = idx( j*4+3, N );
        v3 = idx( j*4+0, N );
        faces.push_back({v0, v1, v2, v3});
        // out << "f " << v0 << " " << v1 << " " << v2 << " " << v3 << endl;

        v0 = idx( i*4+1, N );
        v1 = idx( i*4+2, N );
        v2 = idx( j*4+2, N );
        v3 = idx( j*4+1, N );
        faces.push_back({v0, v1, v2, v3});
        // out << "f " << v0 << " " << v1 << " " << v2 << " " << v3 << endl;

        v0 = idx( j*4+1, N );
        v1 = idx( j*4+3, N );
        v2 = idx( i*4+3, N );
        v3 = idx( i*4+1, N );          
        faces.push_back({v0, v1, v2, v3});
        // out << "f " << v0 << " " << v1 << " " << v2 << " " << v3 << endl;
      }
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
    else if (std::strcmp(poly_str.c_str(), "bunny-hull") == 0){
      std::tie(mesh, geometry) = readManifoldSurfaceMesh("../meshes/bunny-convhull.obj");
      return std::make_tuple(std::move(mesh), std::move(geometry));
    }
    else if (std::strcmp(poly_str.c_str(), "bunnylp") == 0){
      std::tie(mesh, geometry) = readManifoldSurfaceMesh("../meshes/bunny-lp.stl");
      return std::make_tuple(std::move(mesh), std::move(geometry));
    }
    else if (std::strcmp(poly_str.c_str(), "bunnylp-hull") == 0){
      std::tie(mesh, geometry) = readManifoldSurfaceMesh("../meshes/bunny-lp-convhull.obj");
      return std::make_tuple(std::move(mesh), std::move(geometry));
    }
    else if (std::strcmp(poly_str.c_str(), "gomboc") == 0){
      // std::tie(mesh, geometry) = readManifoldSurfaceMesh("../meshes/gomboc2.stl");
      std::tie(mesh, geometry) = readManifoldSurfaceMesh("../meshes/gomboc.obj");
      return std::make_tuple(std::move(mesh), std::move(geometry));
    }
    else if (std::strcmp(poly_str.c_str(), "rndbs 9-gon 1") == 0){
    }
    else if (std::strcmp(poly_str.c_str(), "sphere") == 0){
      std::tie(mesh, geometry) = readManifoldSurfaceMesh("../meshes/sphere.obj");
      return std::make_tuple(std::move(mesh), std::move(geometry));
    }
    else {
      throw std::runtime_error("no valid string provided\n");
    }
    // std::unique_ptr<ManifoldSurfaceMesh> poly_triangulated;
    
    mesh.reset(new ManifoldSurfaceMesh(faces));
    geometry = std::unique_ptr<VertexPositionGeometry>(new VertexPositionGeometry(*mesh));
    for (Vertex v : mesh->vertices()) {
        // Use the low-level indexers here since we're constructing
        (*geometry).inputVertexPositions[v] = positions[v.getIndex()];
    }
    return std::make_tuple(std::move(mesh), std::move(geometry));
}



// generate simple examples
void preprocess_mesh(ManifoldSurfaceMesh* mesh, VertexPositionGeometry* geometry, bool triangulate){
  // readManifoldSurfaceMesh()
  center_and_normalize(mesh, geometry);
  if (triangulate) {
    for (Face f: mesh->faces())
      mesh->triangulate(f);
    mesh->compress();
  }
  for (Edge e: mesh->edges()){
    double dihedangle= geometry->edgeDihedralAngle(e);
    if (dihedangle > PI)
      printf(" edge %d dihedral is large by %f\n", e.getIndex(), dihedangle - PI);
  }
}