#include "arc_algebra.h"

double EPS = 1e-8;
 
Vector3 intersect_arc_ray_with_arc(Vector3 R1, Vector3 R2, Vector3 A, Vector3 B){
    // printf(" norms: %f, %f, %f, %f \n", R1.norm(), R2.norm(), A.norm(), B.norm());
    
    assert(abs(R1.norm() - 1.) <= EPS && abs(R2.norm() - 1.) <= EPS && abs(A.norm() - 1.) <= EPS && abs(B.norm() - 1.) <= EPS);
    Vector3 ray_plane_normal = cross(R1, R2),
            AB_plane_normal = cross(A, B);
    Vector3 intersection_dir = cross(ray_plane_normal, AB_plane_normal);
    Vector3 p1 = intersection_dir.normalize();
    Vector3 p2 = -p1;
    double AB_len = (A-B).norm();
    // can approximate arc len with len for comparison; since arcs are always less than 180 degrees 
    if ((p1 - A).norm() <= AB_len + EPS && (p1 - B).norm() <= AB_len + EPS){ //
        return p1;
    }
    else if ((p2 - A).norm() <= AB_len + EPS && (p2 - B).norm() <= AB_len + EPS){
        return p2;
    }
    else {
        return Vector3::zero();
    }
    
}
