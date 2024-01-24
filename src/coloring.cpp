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
#include "coloring.h"


FaceData<Vector3> generate_random_colors(ManifoldSurfaceMesh *mesh){
  FaceData<Vector3> face_colors(*mesh, Vector3({0.,0.,0.}));
  // printf("mesh size %d, fcolor size: %d\n",mesh->nFaces(), face_colors.size());
  // fingers crossed for a prime enough vector
  Vector3 holy_vec({exp(1.), PI, 1});
  holy_vec = holy_vec.normalize();
  // init point
  Vector3 curr_color = {0.8, 0.6, 0.2}; //{randomReal(0,1), randomReal(0,1), randomReal(0,1)};
  for (Face f: mesh->faces()){
    curr_color += holy_vec;
    if ((curr_color - Vector3::constant(1.)).norm() < 0.1) // too white
      continue;
    curr_color.x = curr_color.x - floor(curr_color.x);
    curr_color.y = curr_color.y - floor(curr_color.y);
    curr_color.z = curr_color.z - floor(curr_color.z);
    face_colors[f] = curr_color;
  }
  return face_colors;
}