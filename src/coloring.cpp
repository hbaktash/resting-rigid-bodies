#include "coloring.h"


FaceData<Vector3> generate_random_colors(ManifoldSurfaceMesh *mesh){
  FaceData<Vector3> face_colors(*mesh);
  // fingers crossed for a prime enough vector
  Vector3 holy_vec({exp(1.), PI, 1});
  holy_vec = holy_vec.normalize();
  // init point
  Vector3 curr_color = {0.8, 0.6, 0.2}; //{randomReal(0,1), randomReal(0,1), randomReal(0,1)};
  for (Face f: mesh->faces()){
    face_colors[f] = curr_color;
    curr_color += holy_vec;
    curr_color.x = curr_color.x - floor(curr_color.x);
    curr_color.y = curr_color.y - floor(curr_color.y);
    curr_color.z = curr_color.z - floor(curr_color.z);
  }
  return face_colors;
}