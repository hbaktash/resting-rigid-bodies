#include "arc_algebra.h"

double EPS = 1e-8;
 

// can approximate arc len with len for comparison; since arcs are always less than 180 degrees 
bool is_on_arc_segment(Vector3 P, Vector3 A, Vector3 B){
    if (P.norm() <= EPS)
        return false;
    double AB_len = (A-B).norm();
    return (P - A).norm() <= AB_len + EPS && (P - B).norm() <= AB_len + EPS;
}

double patch_area(Vector3 A, Vector3 B, Vector3 C, Vector3 D){
    // TODO: double count all traingles and divide by 2 in the end
    return 0;
}

Vector3 intersect_arc_ray_with_arc(Vector3 R1, Vector3 R2, Vector3 A, Vector3 B){
    // printf(" norms: %f, %f, %f, %f \n", R1.norm(), R2.norm(), A.norm(), B.norm());
    
    assert(abs(R1.norm() - 1.) <= EPS && abs(R2.norm() - 1.) <= EPS && abs(A.norm() - 1.) <= EPS && abs(B.norm() - 1.) <= EPS);
    Vector3 ray_plane_normal = cross(R1, R2),
            AB_plane_normal = cross(A, B);
    Vector3 intersection_dir = cross(ray_plane_normal, AB_plane_normal);
    Vector3 p1 = intersection_dir.normalize();
    Vector3 p2 = -p1;
    if (is_on_arc_segment(p1, A, B)) //
        return p1;
    else if (is_on_arc_segment(p2, A, B))
        return p2;
    else 
        return Vector3::zero();
    
}
