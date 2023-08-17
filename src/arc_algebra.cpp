#include "arc_algebra.h"


Vector3 intersect_arc_ray_with_arc(Vector3 R1, Vector3 R2, Vector3 A, Vector3 B){
    assert(R1.norm() == 1. && R2.norm() == 1. && A.norm() == 1. && B.norm() == 1.);
    Vector3 ray_plane_normal = cross(R1, R2);
}
