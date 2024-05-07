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
#include "arc_algebra.h"

double EPS = 1e-8;
 
// can approximate arc len with len for comparison; since arcs are always less than 180 degrees 
bool is_on_arc_segment(Vector3 P, Vector3 A, Vector3 B){
    if (P.norm() <= EPS)
        return false;
    double AB_len = (A-B).norm();
    return (P - A).norm() <= (AB_len + EPS) && (P - B).norm() <= (AB_len + EPS);
}

// portion on A's side
double arc_portion(Vector3 mid_point, Vector3 A, Vector3 B){
    assert(dot(cross(A, B), mid_point) <= EPS); // all on one arc
    double total_angle = angle(A, B),
           PA_angle  = angle(mid_point, A),
           PB_angle  = angle(mid_point, B);
    if (PA_angle >= total_angle && PA_angle >= PB_angle)
        return 1.;
    else if (PB_angle >= total_angle && PB_angle >= PA_angle)
        return 0.;
    else 
        return PA_angle/total_angle;
}

double patch_area(Vector3 A, Vector3 B, Vector3 C, Vector3 D){
    // TODO: double count all traingles and divide by 2 in the end
    // std::cout<< "$$ patch area for \n" <<
    //             "  A:" << A << "\n"<<
    //             "  B:" << B << "\n"<<
    //             "  C:" << C << "\n"<<
    //             "  D:" << D << "\n$$\n";
    double ABC_area = triangle_patch_area_on_sphere(A, B, C),
           ABD_area = triangle_patch_area_on_sphere(A, B, D),
           CDA_area = triangle_patch_area_on_sphere(C, D, A),
           CDB_area = triangle_patch_area_on_sphere(C, D, B);
    // printf("$$ patch areas: %f, %f, %f, %f \n", ABC_area, ABD_area, CDA_area, CDB_area);
    return (ABC_area + ABD_area + CDA_area + CDB_area)/2.;
}

Vector3 intersect_arc_ray_with_arc(Vector3 R1, Vector3 R2, Vector3 A, Vector3 B, bool &sign_change){
    // printf(" norms: %f, %f, %f, %f \n", R1.norm(), R2.norm(), A.norm(), B.norm());
    
    assert(abs(R1.norm() - 1.) <= EPS && abs(R2.norm() - 1.) <= EPS && abs(A.norm() - 1.) <= EPS && abs(B.norm() - 1.) <= EPS);
    Vector3 ray_plane_normal = cross(R1, R2),
            AB_plane_normal = cross(A, B);
    Vector3 intersection_dir = cross(ray_plane_normal, AB_plane_normal);
    Vector3 p1 = intersection_dir.normalize();
    if (is_on_arc_segment(p1, A, B)){
        sign_change = false;
        return p1;
    } 
    else if (is_on_arc_segment(-p1, A, B)){
        sign_change = true;
        return -p1;
    }
    else 
        return Vector3::zero();
    
}



// Gaussian curvature for a general polygonal case
double gaussian_curvature(Vertex v, VertexPositionGeometry &geometry){
    if (v.isBoundary()) return 0.;
    double gaussianCurvature = 2. * PI;
    for (Halfedge he : v.outgoingHalfedges()) {
        Vertex v_next = he.tipVertex(),
               v_pre = he.prevOrbitFace().tailVertex();
        // printf(" Gaussing for pre,v,next: %d, %d, %d\n", v_pre.getIndex(), v.getIndex(), v_next.getIndex());
        Vector3 A = geometry.inputVertexPositions[v],
                B = geometry.inputVertexPositions[v_pre],
                C = geometry.inputVertexPositions[v_next];
        double q = dot(unit(B - A), unit(C - A));
        q = clamp(q, -1.0, 1.0);
        double angle = std::acos(q);
        // printf(" angle : %f\n", angle);
        gaussianCurvature -= angle;
    }
    return gaussianCurvature;
}


// A being in the middle
double angle_on_sphere(Vector3 P1, Vector3 A, Vector3 P2){
    Vector3 n1 = cross(P1, A),
            n2 = cross(A, P2);
    return angle(n1, n2);
}


// area of a triangular patch on sphere
double triangle_patch_area_on_sphere(Vector3 A, Vector3 B, Vector3 C){
    Vector3 n_AB = cross(A, B),
            n_BC = cross(B, C),
            n_CA = cross(C, A);
    if (n_AB.norm() <= EPS || n_BC.norm() <= EPS || n_CA.norm() <= EPS)
        return 0.;
    double angle_A = angle(n_AB, -n_CA),
           angle_B = angle(n_BC, -n_AB),
           angle_C = angle(n_CA, -n_BC);
    return angle_A + angle_B + angle_C - PI;
}


Vector3 project_on_plane(Vector3 point, Vector3 offset, Vector3 normal){
    Vector3 unit_normal = normal.normalize();
    return point + unit_normal * dot(offset - point, unit_normal);
}


Vector3 point_to_segment_normal(Vector3 P, Vector3 A, Vector3 B){
    Vector3 PB = B - P,
            AB = B - A;
    Vector3 ortho_p = PB - AB*dot(AB, PB)/dot(AB,AB);
    return ortho_p;
}




// autodiff version

// autodiff::Vector3var point_to_segment_normal_ad(autodiff::MatrixX3var &poses, autodiff::Vector3var &G, Edge e){
//     Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
//     autodiff::Vector3var PB = poses.row(v2.getIndex()) - G.transpose(),
//                          AB = poses.row(v2.getIndex()) - poses.row(v1.getIndex());
//     return PB - AB * AB.dot(PB)/AB.dot(AB);
// }


// ** only called when intersection has happens
// autodiff::Vector3var intersect_arc_ray_with_arc_ad(autodiff::MatrixX3var &poses, autodiff::Vector3var &G, Vertex v,
//                                                    autodiff::Vector3var &R2, autodiff::Vector3var &A, autodiff::Vector3var &B,
//                                                    bool sign_change){
//     // assertions have been made in earlier non-AD call of intersect_arc_ray_with_arc(...)
//     // i dont know how intermeditate variables are handled in autodiff; so keeping theese function as direct as possible
//     autodiff::Vector3var p = ((poses.row(v.getIndex()) - G.transpose()).cross(R2)).cross(A.cross(B));
//     if (sign_change)
//         return -p;
//     return p;
// }


// autodiff::Vector3var face_normal_ad(autodiff::MatrixX3var &poses, Face f){
//     // Gather vertex positions for next three vertices
//     Halfedge he = f.halfedge();
//     Vertex v1 = he.vertex(),
//            v2 = he.next().vertex(),
//            v3 = he.next().next().vertex();
//     return ((poses.row(v2.getIndex()) - poses.row(v1.getIndex())).cross(poses.row(v3.getIndex()) - poses.row(v1.getIndex())));//.normalized();
// }


autodiff::var triangle_patch_area_on_sphere_ad(autodiff::Vector3var &A, autodiff::Vector3var &B, autodiff::Vector3var &C){
    autodiff::Vector3var n_AB = A.cross(B),
                         n_BC = B.cross(C),
                         n_CA = C.cross(A);
    // if (n_AB.norm() <= EPS || n_BC.norm() <= EPS || n_CA.norm() <= EPS)
    //     return 0.;
    // TODO: use the  det forumla
    // double angle_A = angle(n_AB, -n_CA),
    //        angle_B = angle(n_BC, -n_AB),
    //        angle_C = angle(n_CA, -n_BC);
    // autodiff::var angle_A = acos(n_AB.dot(-n_CA)/(n_AB.norm()*n_CA.norm())),
    //               angle_B = acos(n_BC.dot(-n_AB)/(n_BC.norm()*n_AB.norm())),
    //               angle_C = acos(n_CA.dot(-n_BC)/(n_CA.norm()*n_BC.norm()));
    autodiff::var angle_A = atan2(n_AB.cross(-n_CA).norm(), n_AB.dot(-n_CA)),
                  angle_B = atan2(n_BC.cross(-n_AB).norm(), n_BC.dot(-n_AB)),
                  angle_C = atan2(n_CA.cross(-n_BC).norm(), n_CA.dot(-n_BC));
    return angle_A + angle_B + angle_C - PI;
}


// template <typename DerivedV>
// Eigen::MatrixBase<DerivedV> point_to_segment_normal_ad(Eigen::MatrixBase<DerivedV> &poses, Eigen::VectorX<DerivedV> &P, Edge e){
//     Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
//     Eigen::VectorX<DerivedV> PB = poses[v2.getIndex()] - P,
//                              AB = poses[v2.getIndex()] - poses[v1.getIndex()];
//     return PB - AB * AB.dot(PB)/AB.dot(AB);
// }

// template <typename DerivedV>
// Eigen::MatrixBase<DerivedV> face_normal_autodiff(Eigen::MatrixBase<DerivedV> &poses, Face f){
//     // Gather vertex positions for next three vertices
//     Halfedge he = f.halfedge();
//     Vertex v1 = he.vertex(),
//            v2 = he.next().vertex(),
//            v3 = he.next().next().vertex();
//     return (poses.row(v2.getIndex()) - poses.row(v1.getIndex())).cross(poses.row(v3.getIndex()) - poses.row(v1.getIndex()));
// }

// Derivatives
// std::tuple<Vector3, Vector3, Vector3> point_to_segment_normal_grads(Vector3 P, Vector3 A, Vector3 B){
//     // grad P, grad A, grad B
//     return std::tuple<Vector3, Vector3, Vector3>({Vector3::constant(0.), Vector3::constant(0.), Vector3::constant(0.)});
// }