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

#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"

#include "deformation.h"

#include "geometrycentral/numerical/linear_solvers.h"

// point p with respect to triangle ABC
Vector3 barycentric(Vector3 p, Vector3 A, Vector3 B, Vector3 C) {
    Vector3 v0 = B - A, v1 = C - A, v2 = p - A;
    double d00 = dot(v0, v0),
           d01 = dot(v0, v1),
           d11 = dot(v1, v1),
           d20 = dot(v2, v0),
           d21 = dot(v2, v1);
    double denom = d00 * d11 - d01 * d01;
    Vector3 ans = Vector3::zero();
    ans.y = (d11 * d20 - d01 * d21) / denom;
    ans.z = (d00 * d21 - d01 * d20) / denom;
    ans.x = 1.0f - ans.y - ans.z;
    return ans;
}

// geometrycentral::DenseMatrix<double> get_ARAP_positions(
//                                        geometrycentral::DenseMatrix<double> old_pos_mat,
//                                        geometrycentral::DenseMatrix<double> new_pos_mat, 
//                                        geometrycentral::DenseMatrix<double> init_sol, 
//                                        ManifoldSurfaceMesh &inner_mesh,
//                                        geometrycentral::Vector<int> hull_indices){
//     Eigen::MatrixXd V,U, bc;
//     Eigen::MatrixXi F;
//     Eigen::VectorXi S,b;
//     igl::ARAPData arap_data;
//     arap_data.max_iter = 20;
//     SparseMatrix<double> L;
//     F = inner_mesh.getFaceVertexMatrix<size_t>().cast<int>();
//     V = old_pos_mat;
//     b = hull_indices;
//     // igl::arap_precomputation(old_pos_mat, F, 3, hull_indices, arap_data);
//     igl::arap_precomputation(V, F, 3, b, arap_data);
//     // geometrycentral::DenseMatrix<double> U = init_sol; 
//     U = init_sol;
//     bc = new_pos_mat(hull_indices, Eigen::all);
//     igl::arap_solve(bc, arap_data, U);
//     return U;
// }


DeformationSolver::DeformationSolver(ManifoldSurfaceMesh *_mesh, VertexPositionGeometry *_old_geometry,
                                     ManifoldSurfaceMesh *_convex_mesh, VertexPositionGeometry *_convex_geometry){
    mesh = _mesh;
    old_geometry = _old_geometry;
    convex_mesh = _convex_mesh;
    convex_geometry = _convex_geometry;
}

double DeformationSolver::get_scheduled_weight(double init_w, double final_w){
    double current_val = 0.;
    if (init_w < final_w)
         current_val = init_w/internal_pt;
    else
        current_val = init_w * internal_pt;
    
    if (current_val > std::max(final_w, init_w))
        return std::max(final_w, init_w);
    if (current_val < std::min(final_w, init_w))
        return std::min(final_w, init_w);
    return current_val;
    // return (final_w * init_w) /(init_w + (final_w - init_w) * internal_pt);
}

double DeformationSolver::bending_energy(VertexPositionGeometry *new_geometry){
    assert(new_geometry->mesh.nVertices() == mesh->nVertices()); // should be the same meshes actually

    double energy = 0.;
    for (Edge e: mesh->edges()){
        // dihedral angle difference
        double old_dihedral_angle = old_geometry->edgeDihedralAngle(e),
               new_dihedral_angle = new_geometry->edgeDihedralAngle(e);
        double dihedral_angle_diff = new_dihedral_angle - old_dihedral_angle;
        
        // old geometry quantities
        Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
        Vector3 p1 = old_geometry->inputVertexPositions[v1],
                p2 = old_geometry->inputVertexPositions[v2];
        double e_len_sqrd = (p1 - p2).norm2();
        double area1 = old_geometry->faceArea(e.halfedge().face()),
               area2 = old_geometry->faceArea(e.halfedge().twin().face());
        
        // could ommit the 3. factor
        energy += 3.*dihedral_angle_diff * dihedral_angle_diff * e_len_sqrd / (area1 + area2);
    }
    return energy;
}


VertexData<Vector3> DeformationSolver::bending_energy_gradient(VertexPositionGeometry *new_geometry){
    assert(new_geometry->mesh.nVertices() == mesh->nVertices()); // should be the same meshes actually

    VertexData<Vector3> bending_energy_gradients(*mesh, Vector3::zero());
    

    for (Edge e: mesh->edges()){
        // dihedral angle difference
        double old_dihedral_angle = old_geometry->edgeDihedralAngle(e),
               new_dihedral_angle = new_geometry->edgeDihedralAngle(e);
        double dihedral_angle_diff = new_dihedral_angle - old_dihedral_angle;
        
        // old geometry quantities
        Vertex v1 = e.firstVertex(), v2 = e.secondVertex(),
               u1 = e.halfedge().next().tipVertex(), u2 = e.halfedge().twin().next().tipVertex();
        Vector3 p1 = old_geometry->inputVertexPositions[v1],
                p2 = old_geometry->inputVertexPositions[v2],
                q1 = old_geometry->inputVertexPositions[u1],
                q2 = old_geometry->inputVertexPositions[u2];
        double e_len_sqrd = (p1 - p2).norm2();
        double area1 = old_geometry->faceArea(e.halfedge().face()),
               area2 = old_geometry->faceArea(e.halfedge().twin().face());
        
        // could ommit the 3. factor
        double constant_factor = 3. * dihedral_angle_diff * e_len_sqrd / (area1 + area2);

        // ** assuming outward normals **
        Vector3 u1_theta_grad = dihedral_angle_grad_B(p1, p2, q2, q1),
                u2_theta_grad = dihedral_angle_grad_C(p1, p2, q2, q1),
                v1_theta_grad = dihedral_angle_grad_G(p1, p2, q2, q1),
                v2_theta_grad = dihedral_angle_grad_A(p1, p2, q2, q1);

        bending_energy_gradients[v1] += constant_factor * v1_theta_grad;
        bending_energy_gradients[v2] += constant_factor * v2_theta_grad;
        bending_energy_gradients[u1] += constant_factor * u1_theta_grad;
        bending_energy_gradients[u2] += constant_factor * u2_theta_grad;
    }
    return bending_energy_gradients;
}


double DeformationSolver::closest_point_energy(VertexPositionGeometry *new_geometry){
    // if (closest_point_assignment.size() == 0)
    //     assign_closest_points(new_geometry);
    // assign_closest_points(new_geometry); 
    double energy = 0.;
    DenseMatrix<double> new_pos_mat = vertex_data_to_matrix(new_geometry->inputVertexPositions),
                        convex_points_mat = vertex_data_to_matrix(convex_geometry->inputVertexPositions);
    DenseMatrix<double> closest_assignments = closest_point_operator*new_pos_mat;
    energy = (closest_assignments - convex_points_mat).norm();
    // polyscope::registerPointCloud("closest assignments", closest_assignments);
    // polyscope::show();
    // for (Vertex c_v: convex_mesh->vertices()){
    //     Vector3 c_p = convex_geometry->inputVertexPositions[c_v],
    //             closest_point = closest_point_assignment[c_v].interpolate(new_geometry->inputVertexPositions);
    //     energy += (c_p - closest_point).norm();
    // }
    // printf(" CP energy: %f \n", energy);
    return energy;
}

double DeformationSolver::closest_point_energy(Vector<double> flat_new_pos_mat){
    DenseMatrix<double> convex_points_mat = vertex_data_to_matrix(convex_geometry->inputVertexPositions);
    return (closest_point_flat_operator*flat_new_pos_mat - tinyAD_flatten(convex_points_mat)).norm();
}

DenseMatrix<double> DeformationSolver::closest_point_energy_gradient(VertexPositionGeometry *new_geometry){
    VertexData<Vector3> CP_energy_gradients(*mesh, Vector3::zero());
    DenseMatrix<double> new_pos_mat = vertex_data_to_matrix(new_geometry->inputVertexPositions),
                        convex_points_mat = vertex_data_to_matrix(convex_geometry->inputVertexPositions);
    DenseMatrix<double> CP_grad = closest_point_operator.transpose()*
                                  (closest_point_operator*new_pos_mat - convex_points_mat);
    return CP_grad;
}


void DeformationSolver::assign_closest_vertices(VertexPositionGeometry *new_geometry, bool allow_multi_assignment){
    closest_point_assignment = VertexData<SurfacePoint>(*convex_mesh, SurfacePoint());
    closest_point_distance = VertexData<double>(*convex_mesh, -1.);
    CP_involvement = VertexData<bool>(*mesh, false);
    // stats
    size_t vertex_double_hits = 0;
    std::vector<Vector3> double_hit_CVs, double_hit_IPs;

    // linear operator corresponding to assignments
    closest_point_operator = Eigen::SparseMatrix<double>(convex_mesh->nVertices(), mesh->nVertices());
    closest_point_flat_operator = Eigen::SparseMatrix<double>(3 * convex_mesh->nVertices(), 3*mesh->nVertices());
    std::vector<Eigen::Triplet<double>> tripletList, flat_tripletList;
    tripletList.reserve(3*mesh->nVertices());
    flat_tripletList.reserve(9*mesh->nVertices());

    for (Vertex c_v: convex_mesh->vertices()){
        Vector3 c_p = convex_geometry->inputVertexPositions[c_v];

        Vertex closest_vertex;
        double min_vertex_dist  = std::numeric_limits<double>::infinity();
        bool last_hit_was_taken = false;
        Vector3 last_taken_point;
        for (Vertex v: mesh->vertices()){
            Vector3 A = new_geometry->inputVertexPositions[v];
            double vertex_dist = (A - c_p).norm();
            if (vertex_dist < min_vertex_dist){
                if (CP_involvement[v]){
                    last_hit_was_taken = true;
                    last_taken_point = A;
                    if (!allow_multi_assignment) continue; // go for the next closest point
                    // otherwise a point can be assigned as multiple closest points 
                }
                closest_vertex = v;
                min_vertex_dist = vertex_dist;
                if (!allow_multi_assignment) last_hit_was_taken = false;
            }            
        }
        if (last_hit_was_taken){
            vertex_double_hits++;
            // DEBUG
            double_hit_CVs.push_back(c_p);
            double_hit_IPs.push_back(last_taken_point);
        }
        
        // update assignment, involvement, distance
        closest_point_assignment[c_v] = SurfacePoint(closest_vertex); 
        closest_point_distance[c_v] = min_vertex_dist;
        CP_involvement[closest_vertex] = true;

        // sparse matrix insertion
        tripletList.emplace_back(c_v.getIndex(), closest_vertex.getIndex(), 1.);
        for (size_t d = 0; d < 3; d++)
            flat_tripletList.emplace_back(3 * c_v.getIndex() + d, 3*closest_vertex.getIndex() + d, 1.);
    }
    closest_point_operator.setFromTriplets(tripletList.begin(), tripletList.end());
    closest_point_flat_operator.setFromTriplets(flat_tripletList.begin(), flat_tripletList.end());
    // auto dp_cvs = polyscope::registerPointCloud(" double hit CVs ", double_hit_CVs);
    // auto dp_IPs = polyscope::registerPointCloud(" double hit interiors ", double_hit_IPs);
    // dp_cvs->setPointColor({0.8,0.3,0.3});
    // dp_IPs->setPointColor({0.1,0.1,0.9});
    printf("  ** stats for vertex-only CP: double hits %d \n", vertex_double_hits);
}

void DeformationSolver::assign_closest_points_barycentric(VertexPositionGeometry *new_geometry){
    threshold_exceeded = 0;
    marked_to_split = VertexData<bool>(*convex_mesh, false);
    closest_point_assignment = VertexData<SurfacePoint>(*convex_mesh, SurfacePoint());
    closest_point_distance = VertexData<double>(*convex_mesh, -1.);
    CP_involvement = VertexData<bool>(*mesh, false);
    // stats
    face_assignments = 0; edge_assignments = 0; vertex_assignments = 0;

    // linear operator corresponding to assignments
    closest_point_operator = Eigen::SparseMatrix<double>(convex_mesh->nVertices(), mesh->nVertices());
    closest_point_flat_operator = Eigen::SparseMatrix<double>(3 * convex_mesh->nVertices(), 3*mesh->nVertices());
    std::vector<Eigen::Triplet<double>> tripletList, flat_tripletList;
    tripletList.reserve(3*mesh->nVertices());
    flat_tripletList.reserve(9*mesh->nVertices());
    // DEBUG
    // std::vector<Vector3> barry_points, bary_normals;
    // for (Face f: mesh->faces()){
    //   Halfedge he = f.halfedge();
    //   Vector3 barry = (new_geometry->inputVertexPositions[he.vertex()] + new_geometry->inputVertexPositions[he.next().vertex()] + new_geometry->inputVertexPositions[he.next().next().vertex()])/3.;
    //   barry_points.push_back(barry);
    //   bary_normals.push_back(new_geometry->faceNormal(f));
    // }
    // auto bary_pc = polyscope::registerPointCloud("barycenters", barry_points);
    // bary_pc->addVectorQuantity("normals", bary_normals);
    // polyscope::show();
    for (Vertex c_v: convex_mesh->vertices()){
        if (frozen_assignment[c_v].getIndex() != INVALID_IND)
            continue;
        Face closest_face;
        Edge closest_edge;
        Vertex closest_vertex;
        
        Vector3 on_edge_projection;
        double min_face_dist  = std::numeric_limits<double>::infinity(),
               min_edge_dist   = std::numeric_limits<double>::infinity(),
               min_vertex_dist  = std::numeric_limits<double>::infinity();

        Vector3 c_p = convex_geometry->inputVertexPositions[c_v];

        if (!vertex_only_assignment[c_v]){
            // there should be a smarter way of checking all elements in one loop; maybe even without a flagging them??
            for (Face f: mesh->faces()){ 
                Vector3 f_normal = new_geometry->faceNormal(f);
                Halfedge curr_he  = f.halfedge(),
                        first_he = f.halfedge();
                // assume outward normals??
                double face_dist = dot(f_normal, c_p - new_geometry->inputVertexPositions[f.halfedge().vertex()]); // using some point of f
                if (face_dist < 0) continue; // not on the right side of the face (outward normal)
                bool face_is_projectable = true;
                while (true){ // checking if face-projectable
                    Vertex v1 = curr_he.tailVertex(), v2 = curr_he.tipVertex();
                    Vector3 A = new_geometry->inputVertexPositions[v1], 
                            B = new_geometry->inputVertexPositions[v2];
                    Vector3 N_PAB = cross(A - c_p, B - A);
                    if (dot(f_normal, N_PAB) <= 0) { // not face-projectable on this face
                        face_is_projectable = false;
                        break; // go to next face
                    }
                    curr_he = curr_he.next();
                    if (curr_he == first_he)
                        break;
                }
                // if (face_is_projectable) polyscope::warning("some face is projectable");
                if (face_is_projectable && face_dist < min_face_dist){
                    min_face_dist = face_dist;
                    closest_face = f;
                }
            }
            for (Edge e: mesh->edges()){
                Vector3 A = new_geometry->inputVertexPositions[e.firstVertex()], 
                        B = new_geometry->inputVertexPositions[e.secondVertex()];
                Vector3 PB = B - c_p,
                        PA = A - c_p,
                        AB = B - A;
                bool edge_is_projectable = dot(PB, AB) >= 0 && dot(PA, -AB) >= 0;
                if (edge_is_projectable){
                    Vector3 ortho_p = PB - AB*dot(AB, PB)/dot(AB,AB);
                    double edge_dist =  ortho_p.norm();
                    if (edge_dist < min_edge_dist){
                        closest_edge = e;
                        min_edge_dist = edge_dist;
                        on_edge_projection = c_p + ortho_p;
                    }
                }
            }
        }
        for (Vertex v: mesh->vertices()){
            Vector3 A = new_geometry->inputVertexPositions[v];
            double vertex_dist = (A - c_p).norm();
            if (vertex_dist < min_vertex_dist){
                closest_vertex = v;
                min_vertex_dist = vertex_dist;
            }
        }
        

        // assign SurfacePoint and assign barycentric coordinates
        // printf(" at cv %d fd %f, ed %f, vd %f\n", c_v.getIndex(), min_face_dist, min_edge_dist, min_vertex_dist);
        double min_dist = std::min(min_vertex_dist, std::min(min_edge_dist, min_face_dist)); // ordering is important in case of equality
        
        // SurfacePoint assignment
        if (min_dist == min_face_dist){ // assuming triangle mesh
            Vector3 f_normal = new_geometry->faceNormal(closest_face);
            Vector3 on_face_projection = c_p - f_normal * dot(f_normal, 
                                                              c_p - new_geometry->inputVertexPositions[closest_face.halfedge().vertex()]);
            Vertex v1 = closest_face.halfedge().vertex(),
                   v2 = closest_face.halfedge().next().vertex(),
                   v3 = closest_face.halfedge().next().next().vertex();
            Vector3 A = new_geometry->inputVertexPositions[v1],
                    B = new_geometry->inputVertexPositions[v2],
                    C = new_geometry->inputVertexPositions[v3];
            Vector3 bary_coor = barycentric(on_face_projection, A, B, C);
            closest_point_assignment[c_v] = get_robust_barycentric_point(SurfacePoint(closest_face, bary_coor), split_robustness_threshold);
        }
        else if (min_dist == min_edge_dist){
            Vertex v1 = closest_edge.firstVertex(),
                   v2 = closest_edge.secondVertex();
            Vector3 A = new_geometry->inputVertexPositions[v1], 
                    B = new_geometry->inputVertexPositions[v2];
            double tVal = (on_edge_projection - A).norm()/(A - B).norm();
            closest_point_assignment[c_v] = get_robust_barycentric_point(SurfacePoint(closest_edge, tVal), split_robustness_threshold); 
        }
        else {// min_dist == min_vertex_dist
            if (vertex_only_assignment[c_v] && min_dist <= refinement_CP_threshold && enforce_snapping){
                frozen_assignment[c_v] = closest_vertex;
                closest_point_distance[c_v] = 0.;
                continue;
            }
            closest_point_assignment[c_v] = SurfacePoint(closest_vertex); 
        }

        // Matrix filling and splitting decision
        SurfacePoint robustSP = closest_point_assignment[c_v];
        if (robustSP.type == SurfacePointType::Face){ // FACE
            face_assignments++;
            Face f = robustSP.face;
            Vertex v1 = f.halfedge().vertex(),
                   v2 = f.halfedge().next().vertex(),
                   v3 = f.halfedge().next().next().vertex();
            Vector3 A = new_geometry->inputVertexPositions[v1],
                    B = new_geometry->inputVertexPositions[v2],
                    C = new_geometry->inputVertexPositions[v3];
            double d1 = get_point_distance_to_convex_hull(vec32vec(A)), 
                   d2 = get_point_distance_to_convex_hull(vec32vec(B)), 
                   d3 = get_point_distance_to_convex_hull(vec32vec(C));
            if (d1 <= refinement_CP_threshold && d2 <= refinement_CP_threshold && d3 <= refinement_CP_threshold) {// or max?
                marked_to_split[c_v] = true;
                threshold_exceeded++;
            }
            // update CP involvement
            CP_involvement[v1] = true;
            CP_involvement[v2] = true;
            CP_involvement[v3] = true;
            // operator entries
            Vector3 bary_coor = robustSP.faceCoords;
            tripletList.emplace_back(c_v.getIndex(), v1.getIndex(), bary_coor.x);
            tripletList.emplace_back(c_v.getIndex(), v2.getIndex(), bary_coor.y);
            tripletList.emplace_back(c_v.getIndex(), v3.getIndex(), bary_coor.z);
            for (size_t d = 0; d < 3; d++){
                flat_tripletList.emplace_back(3 * c_v.getIndex() + d, 3*v1.getIndex() + d, bary_coor.x);
                flat_tripletList.emplace_back(3 * c_v.getIndex() + d, 3*v2.getIndex() + d, bary_coor.y);
                flat_tripletList.emplace_back(3 * c_v.getIndex() + d, 3*v3.getIndex() + d, bary_coor.z);
            }
        }
        else if (robustSP.type == SurfacePointType::Edge){ // EDGE
            edge_assignments++;
            Edge e = robustSP.edge;
            Vertex v1 = e.firstVertex(),
                   v2 = e.secondVertex();
            Vector3 A = new_geometry->inputVertexPositions[v1], 
                    B = new_geometry->inputVertexPositions[v2];
            double d1 = get_point_distance_to_convex_hull(vec32vec(A)), 
                    d2 = get_point_distance_to_convex_hull(vec32vec(B));
            if (d1 <= refinement_CP_threshold && d2 <= refinement_CP_threshold){ // close and not robust
                marked_to_split[c_v] = true;
                threshold_exceeded++;
            }
            // update CP involvement
            CP_involvement[v1] = true;
            CP_involvement[v2] = true;
            // operator entries
            double tVal = robustSP.tEdge;
            tripletList.emplace_back(c_v.getIndex(), v1.getIndex(), 1. - tVal);
            tripletList.emplace_back(c_v.getIndex(), v2.getIndex(), tVal);
            for (size_t d = 0; d < 3; d++){
                flat_tripletList.emplace_back(3 * c_v.getIndex() + d, 3*v1.getIndex() + d, 1. - tVal);
                flat_tripletList.emplace_back(3 * c_v.getIndex() + d, 3*v2.getIndex() + d, tVal);
            }
        }
        if (robustSP.type == SurfacePointType::Vertex){//VERTEX
            vertex_assignments++;
            Vertex v = robustSP.vertex;
            // update CP involvement
            CP_involvement[v] = true;
            tripletList.emplace_back(c_v.getIndex(), v.getIndex(), 1.);
            for (size_t d = 0; d < 3; d++)
                flat_tripletList.emplace_back(3 * c_v.getIndex() + d, 3*v.getIndex() + d, 1.);
        }

        closest_point_distance[c_v] = min_dist;
    }
    closest_point_operator.setFromTriplets(tripletList.begin(), tripletList.end());
    closest_point_flat_operator.setFromTriplets(flat_tripletList.begin(), flat_tripletList.end());
    printf("  ** stats for CP: vertices %d, edges: %d, faces: %d \n", vertex_assignments, edge_assignments, face_assignments);
}


Eigen::VectorX<bool> DeformationSolver::get_frozen_flags() {
    Eigen::VectorX<bool> frozen_flags(mesh->nVertices());
    frozen_flags.setConstant(false);
    for (Vertex v: convex_mesh->vertices()){
        if (frozen_assignment[v].getIndex() != INVALID_IND){
            size_t i = frozen_assignment[v].getIndex();
            frozen_flags(3*i + 0) = true;
            frozen_flags(3*i + 1) = true;
            frozen_flags(3*i + 2) = true;
        }
    }
    return frozen_flags;
}


Eigen::VectorXd DeformationSolver::get_frozen_x(){
    Eigen::VectorXd frozen_x = Eigen::VectorXd::Zero(3*mesh->nVertices());
    for (Vertex cv: convex_mesh->vertices()){
        Vertex inner_assignment = frozen_assignment[cv];
        if (inner_assignment.getIndex() != INVALID_IND){
            size_t i = inner_assignment.getIndex();
            Vector3 cp = convex_geometry->inputVertexPositions[cv];
            frozen_x[3*i + 0] = cp.x;
            frozen_x[3*i + 1] = cp.y;
            frozen_x[3*i + 2] = cp.z;
        }
    }
    return frozen_x;
}


SurfacePoint DeformationSolver::get_robust_barycentric_point(SurfacePoint p, double threshold){
    if (p.type == SurfacePointType::Face){
        Vector3 bary_coords = p.faceCoords;
        bool remove_x = bary_coords.x < threshold,
             remove_y = bary_coords.y < threshold,
             remove_z = bary_coords.z < threshold;
        if      (remove_x && !remove_y && !remove_z)
            return SurfacePoint(p.face.halfedge().next(), 1. - bary_coords.y); 
        else if (!remove_x && remove_y && !remove_z)
            return SurfacePoint(p.face.halfedge().next().next(), 1. - bary_coords.z);
        else if (!remove_x && !remove_y && remove_z)
            return SurfacePoint(p.face.halfedge(), 1. - bary_coords.x);
        else if (remove_x && remove_y && !remove_z)
            return SurfacePoint(p.face.halfedge().next().tipVertex()); // vertex 2
        else if (!remove_x && remove_y && remove_z)
            return SurfacePoint(p.face.halfedge().vertex()); // vertex 0
        else if (remove_x && !remove_y && remove_z)
            return SurfacePoint(p.face.halfedge().next().tailVertex()); // vertex 1
        else if (!remove_x && !remove_y && !remove_z) // all-remove not possible
            return p; 
        else throw std::logic_error(" barycenteric weights should sum to one!");
    }
    else if (p.type == SurfacePointType::Edge){
        if (p.tEdge < threshold)
            return SurfacePoint(p.edge.firstVertex());
        else if (1. - p.tEdge < threshold)
            return SurfacePoint(p.edge.secondVertex());
        else 
            return p;
    }
    else {
        return p;
    }
}


bool DeformationSolver::split_barycentric_closest_points(VertexPositionGeometry *new_geometry){  
    // assign CP should be called already; containers initialized
    size_t mutli_edge_hits = 0, multi_face_hits = 0, multi_vertex_hits = 0;
    EdgeData<bool> edge_is_hit(*mesh, false);
    FaceData<bool> face_is_hit(*mesh, false);
    VertexData<bool> vertex_is_hit(*mesh, false);
    bool split_occured = false;
    for (Vertex cv: convex_mesh->vertices()){
        if (!marked_to_split[cv])
            continue;
        // if (closest_point_distance[cv] > refinement_CP_threshold && local_split)
        //     continue;
        vertex_only_assignment[cv] = true; // TODO make this an option
        split_occured = true;
        SurfacePoint assigned_cp = closest_point_assignment[cv];
        SurfacePoint robust_sp = get_robust_barycentric_point(assigned_cp, split_robustness_threshold); // hitting twice with robustness; could go face->vertex
        if (robust_sp.type != assigned_cp.type)
            polyscope::registerPointCloud("robusted SP", std::vector<Vector3>({assigned_cp.interpolate(new_geometry->inputVertexPositions) ,robust_sp.interpolate(new_geometry->inputVertexPositions)}))->setPointColor({0.8,0.1,0.1});
        Vertex new_v;
        Vector3 split_p_old_geo,
                split_p_new_geo;

        // TODO: robust barycenter check
        if (robust_sp.type == SurfacePointType::Face){
            if (!face_is_hit[robust_sp.face])
                face_is_hit[robust_sp.face] = true;
            else {
                multi_face_hits++;
                vertex_only_assignment[cv] = false;
                continue;
            } 
            // printf("splitting face %d\n", robust_sp.face.getIndex());
            // must be done before spliting
            split_p_old_geo = robust_sp.interpolate(old_geometry->inputVertexPositions);
            split_p_new_geo = robust_sp.interpolate(new_geometry->inputVertexPositions);
            
            new_v = mesh->insertVertex(robust_sp.face);
        }
        else if (robust_sp.type == SurfacePointType::Edge){
            if (!edge_is_hit[robust_sp.edge])
                edge_is_hit[robust_sp.edge] = true;
            else { // not handling multi edge hits at the moment
                mutli_edge_hits++;
                vertex_only_assignment[cv] = false;
                continue;
            }
            // printf("splitting edge %d: %d,%d \n", robust_sp.edge.getIndex(), robust_sp.edge.firstVertex().getIndex(), robust_sp.edge.secondVertex().getIndex());
            // must be done before spliting
            split_p_old_geo = robust_sp.interpolate(old_geometry->inputVertexPositions);
            split_p_new_geo = robust_sp.interpolate(new_geometry->inputVertexPositions);
            // debug 
            new_v = mesh->splitEdgeTriangular(robust_sp.edge).vertex();
        }
        else {
            // printf("splitting vertex %d\n", robust_sp.vertex.getIndex());
            // must be done before spliting
            split_p_old_geo = robust_sp.interpolate(old_geometry->inputVertexPositions);
            split_p_new_geo = robust_sp.interpolate(new_geometry->inputVertexPositions);
            new_v = robust_sp.vertex;
            if (!vertex_is_hit[robust_sp.vertex])
                vertex_is_hit[robust_sp.vertex] = true;
            else 
                multi_vertex_hits++;
        }
        old_geometry->inputVertexPositions[new_v] = split_p_old_geo;
        new_geometry->inputVertexPositions[new_v] = split_p_new_geo;
    }
    mesh->compress();
    if (mutli_edge_hits + multi_face_hits + multi_vertex_hits > 0)
        printf(" splitting is over. multi edge hits %d, \n\t\t\t\t multi face hits: %d, \n\t\t\t\t multi vertex hits: %d\n", mutli_edge_hits, multi_face_hits, multi_vertex_hits);
    return split_occured;
}


Vector<double> DeformationSolver::get_CP_barrier_multiplier_vetor(double CP_involed_constant){
    Vector<double> barrier_multipliers = Vector<double>::Ones(mesh->nVertices());
    for (Vertex v: mesh->vertices()){
        if (CP_involvement[v]){
            barrier_multipliers[v.getIndex()] = CP_involed_constant;
        }
    }
    return barrier_multipliers;
}


void DeformationSolver::build_constraint_matrix_and_rhs(){
    size_t nf = convex_mesh->nFaces();
    constraint_matrix = DenseMatrix<double>::Zero(nf, 3);
    constraint_rhs = Vector<double>::Zero(nf);
    for (Face f: convex_mesh->faces()){
        Vector3 face_normal = convex_geometry->faceNormal(f);
        constraint_matrix.row(f.getIndex()) = vec32vec(face_normal);
        Vector3 point_on_face_normal = convex_geometry->inputVertexPositions[f.halfedge().vertex()];
        constraint_rhs(f.getIndex()) = dot(face_normal, point_on_face_normal);
    }
}


std::tuple<double, DenseMatrix<double>, std::vector<DenseMatrix<double>>> 
DeformationSolver::get_log_barrier_stuff(DenseMatrix<double> new_pos_mat, Vector<double> weights){
    double energy = 0.;
    size_t n = new_pos_mat.rows(),
           nf = convex_mesh->nFaces();
    assert(new_pos_mat.cols() == 3);
    DenseMatrix<double> NP = constraint_matrix * new_pos_mat.transpose(); // nf by n; col j is N*P_j
    DenseMatrix<double> repeated_rhs = constraint_rhs.replicate(1, n);
    assert(repeated_rhs.cols() == n); // same as the number of new pos vertices
    DenseMatrix<double> diff_ij = repeated_rhs - NP; // col j is rhs - N*P_j; i.e. -f_i(P_j) which should be positive
    // slow?
    assert((diff_ij.array() > DenseMatrix<double>::Zero(diff_ij.rows(), diff_ij.cols()).array()).all()); // check interior-ness before entering this function?
    // slow?
    DenseMatrix<double> energy_ij = diff_ij.array().log().matrix(); // check the cost for this; why should I .array() it??
    energy = -energy_ij.sum();
    
    // gradient
    DenseMatrix<double> grads = DenseMatrix<double>::Zero(n, 3);
    for (size_t j = 0; j < n; j++){
        grads.row(j) = (diff_ij.col(j).cwiseInverse().asDiagonal() * constraint_matrix).colwise().sum(); // sum over N_i / (diff_ij)
    }

    // Hessian is a diagonal matrix; so using matrix per vertex
    std::vector<DenseMatrix<double>> hessians(n);
    for (size_t j = 0; j < n; j++){
        hessians[j] = constraint_matrix.transpose() * diff_ij.col(j).cwiseAbs2().cwiseInverse().asDiagonal() * constraint_matrix;
        // hessians[j] = DenseMatrix<double>::Zero(3,3);
        // for (size_t i = 0; i < nf; i++){
        //     DenseMatrix<double> g_gT = grads.row(j).transpose() * grads.row(j);
        //     if(!(g_gT.cols() == 3 && g_gT.rows() == 3))
        //         throw std::logic_error(" gt size "+ std::to_string(g_gT.rows()) + "," + std::to_string(g_gT.cols()) + "\n");
        //     double diff_ij_sqr = diff_ij.coeff(i,j) * diff_ij.coeff(i,j);
        //     hessians[j] += g_gT/diff_ij_sqr; // sum over N_i / (diff_ij)
        // }
    }

    return std::tuple<double, DenseMatrix<double>, std::vector<DenseMatrix<double>>>(energy, grads, hessians);
}


double DeformationSolver::get_log_barrier_energy(DenseMatrix<double> new_pos_mat){
    size_t n = new_pos_mat.rows();
    DenseMatrix<double> NP = constraint_matrix * new_pos_mat.transpose(); // nf by n; col j is N*P_j
    DenseMatrix<double> repeated_rhs = constraint_rhs.replicate(1, n);
    assert(repeated_rhs.cols() == n); // same as the number of new pos vertices
    DenseMatrix<double> diff_ij = repeated_rhs - NP; // col j is N*P_j - rhs; i.e. -f_i(P_j) which should be positive
    DenseMatrix<double> energy_ij = diff_ij.array().log().matrix(); // check the cost for this; why should I .array() it??

    return -energy_ij.sum();
}


bool DeformationSolver::check_feasibility(DenseMatrix<double> new_pos_mat){
    size_t n = new_pos_mat.rows();
    DenseMatrix<double> NP = constraint_matrix * new_pos_mat.transpose(); // nf by n; col j is N*P_j
    DenseMatrix<double> repeated_rhs = constraint_rhs.replicate(1, n);
    assert(repeated_rhs.cols() == n); // same as the number of new pos vertices
    DenseMatrix<double> diff_ij = repeated_rhs - NP; // col j is N*P_j - rhs; i.e. -f_i(P_j) which should be positive
    return (diff_ij.array() > DenseMatrix<double>::Zero(diff_ij.rows(), diff_ij.cols()).array()).all();
}


DenseMatrix<bool> DeformationSolver::get_active_set_matrix(DenseMatrix<double> new_pos_mat, double active_threshold){
    size_t n = new_pos_mat.rows();
    DenseMatrix<double> NP = constraint_matrix * new_pos_mat.transpose(); // nf by n; col j is N*P_j
    DenseMatrix<double> repeated_rhs = constraint_rhs.replicate(1, n);
    assert(repeated_rhs.cols() == n); // same as the number of new pos vertices
    DenseMatrix<double> diff_ij = repeated_rhs - NP; // col j is N*P_j - rhs; all positive
    DenseMatrix<bool> active_constraints = diff_ij.array() < DenseMatrix<double>::Constant(constraint_matrix.rows(), n, active_threshold).array();
    return active_constraints;
}


double DeformationSolver::get_point_distance_to_convex_hull(Eigen::VectorXd p, bool from_faces){
    double distance = 0.;
    if (from_faces){
        size_t m = convex_mesh->nFaces();
        Eigen::VectorXd Np = constraint_matrix * p;
        assert(Np.size() == m); // p or pT?
        Eigen::MatrixXd diff_ij = constraint_rhs - Np; // col j is N*P_j - rhs; i.e. -f_i(P_j) which should be positive
        distance = diff_ij.minCoeff();
    }
    else { // from vertices
        size_t nc = convex_mesh->nVertices(); 
        Eigen::MatrixXd repeated_p = p.transpose().replicate(nc, 1); // DONT FORGET THE TRANSPOSE EVER AGAIN
        assert(repeated_p.rows() == nc && repeated_p.cols() == 3);
        Eigen::MatrixXd diff_vecs = repeated_p - vertex_data_to_matrix(convex_geometry->inputVertexPositions);
        Eigen::VectorXd norms = diff_vecs.rowwise().norm();
        distance = norms.minCoeff();
    }
    return distance;
}


void DeformationSolver::update_bending_rest_constants(){
    // rest dihedrals
    old_geometry->requireEdgeDihedralAngles();
    rest_dihedral_angles = old_geometry->edgeDihedralAngles;
    old_geometry->unrequireEdgeDihedralAngles();

    // rest e/h_e
    rest_bending_constant = EdgeData<double>(*mesh, 0.); // e/h_e
    old_geometry->refreshQuantities();
    for (Edge e : mesh->edges()) {
        Vector3 p1 = old_geometry->inputVertexPositions[e.firstVertex()];
        Vector3 p2 = old_geometry->inputVertexPositions[e.secondVertex()];
        // length and area
        double e_len_sqrd = (p1 - p2).norm2();
        double area1 = old_geometry->faceArea(e.halfedge().face()),
               area2 = old_geometry->faceArea(e.halfedge().twin().face());
        rest_bending_constant[e] = e_len_sqrd/(area1+area2);
        if (rest_bending_constant[e] >= 1e9){ // DEBUG
            printf(" edge Death: %d\n", e.isDead());
            printf("rest constant: %d: %d, %d  = %f\n", e.getIndex(), e.firstVertex().getIndex(), e.secondVertex().getIndex(), rest_bending_constant[e]);
            printf(" -- areas: %f, %f\n", area1, area2);
            printf(" -- length sqrd: %f\n", e_len_sqrd);
            Vector3 p3 = old_geometry->inputVertexPositions[e.halfedge().next().tipVertex()],
                    p4 = old_geometry->inputVertexPositions[e.halfedge().twin().next().tipVertex()];
            std::cout << " -- p1: " << p1 << std::endl;
            std::cout << " -- p2: " << p2 << std::endl;
            std::cout << " -- p3: " << p3 << std::endl;
            std::cout << " -- p4: " << p4 << std::endl;
            auto debug_ps = polyscope::registerSurfaceMesh("debug mesh", old_geometry->inputVertexPositions, mesh->getFaceVertexList());
            FaceData<Vector3> fcolors(*mesh, Vector3::constant(1.));
            fcolors[e.halfedge().face()] = Vector3({1.,0,0});
            fcolors[e.halfedge().twin().face()] = Vector3({1.,0,0});
            debug_ps->addFaceColorQuantity("zero faces!?", fcolors);
            auto debug_pc = polyscope::registerPointCloud("zero areas quad", std::vector<Vector3>({p1, p2, p3, p4}));
            polyscope::show();
        }
    };
}


void DeformationSolver::update_membrane_rest_constants(){
    // rest face areas
    old_geometry->requireFaceAreas();
    rest_face_areas = old_geometry->faceAreas;
    
    // rest metric tensor
    rest_membrane_I_inverted = FaceData<Eigen::Matrix2d>(*mesh);
    for (Face f: mesh->faces()){
        Vertex v1 = f.halfedge().tailVertex(),
               v2 = f.halfedge().tipVertex(),
               v3 = f.halfedge().next().tipVertex();
        Vector3 p1 = old_geometry->inputVertexPositions[v1], 
                p2 = old_geometry->inputVertexPositions[v2],
                p3 = old_geometry->inputVertexPositions[v3];
        Vector3 e2 = p2 - p1, 
                e3 = p3 - p1;
        // Eigen::Matri2<double, 3, 2> J = TinyAD::col_mat(p2 - p1, p3 - p1);
        Eigen::Matrix2d I { // J^T J
            {dot(e2,e2), dot(e2,e3)},
            {dot(e3,e2), dot(e3,e3)}
        };
        rest_membrane_I_inverted[f] = I.inverse();
    }
}


auto DeformationSolver::get_tinyAD_bending_function(){
    // bending function
    auto bendingEnergy_func = TinyAD::scalar_function<3>(mesh->vertices()); // 
    // Add objective term per edge
    // double edge_count = 1.;//(double) mesh->nEdges();
    bendingEnergy_func.add_elements<4>(mesh->edges(), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {
        // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);

        // Get variables 3D vertex positions
        Edge e = element.handle;
        
        if (e.isBoundary()) return (T)0.0;

        Eigen::Vector3<T> p1 = element.variables(e.firstVertex());
        Eigen::Vector3<T> p2 = element.variables(e.secondVertex());
        
        Eigen::Vector3<T> p3 = element.variables(e.halfedge().next().tipVertex());
        Eigen::Vector3<T> p4 = element.variables(e.halfedge().twin().next().tipVertex());

        Eigen::Vector3<T> N1 = (p2 - p1).cross(p3 - p1); //faceNormals[e.halfedge().face()];
        Eigen::Vector3<T> N2 = (p1 - p2).cross(p4 - p2);
        Eigen::Vector3<T> edgeDir = (p2 - p1).normalized();
        T dihedral_angle = atan2(edgeDir.dot(N1.cross(N2)), N1.dot(N2));
        // if (dihedral_angle > rest_dihedral_angles[e]) // was trying smth funny here
        //     return (dihedral_angle - rest_dihedral_angles[e]) * sqrt(rest_bending_constant[e]);
        // else 
        //     return (rest_dihedral_angles[e] - dihedral_angle) * sqrt(rest_bending_constant[e]);
        // T current_edge_len = (p2 - p1).norm();
        return (dihedral_angle - rest_dihedral_angles[e]) * (dihedral_angle - rest_dihedral_angles[e]) * rest_bending_constant[e];
    });

    return bendingEnergy_func;
}


auto DeformationSolver::get_tinyAD_membrane_function(){
    double membrane_mu = 0.2, membrane_lambda = 0.1;
    auto membraneEnergy_func = TinyAD::scalar_function<3>(mesh->vertices()); // 
    membraneEnergy_func.add_elements<3>(mesh->faces(), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {
        // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);

        // Get variables 3D vertex positions
        Face f = element.handle;

        Eigen::Vector3<T> p1 = element.variables(f.halfedge().tailVertex());
        Eigen::Vector3<T> p2 = element.variables(f.halfedge().tipVertex());
        Eigen::Vector3<T> p3 = element.variables(f.halfedge().next().tipVertex());
        Eigen::Matrix<T, 3, 2> J = TinyAD::col_mat(p2 - p1, p3 - p1);
        Eigen::Matrix2<T> I = J.transpose() * J;
        
        Eigen::Matrix2<T> simil_matrix = rest_membrane_I_inverted[f] * I;
        // W(I^{-1} I_tilde)
        
        // T current_triangle_area = 0.5 * (p2 - p1).cross(p3 - p1).norm();
        return rest_face_areas[f] * ((simil_matrix.transpose() * simil_matrix).trace()*0.5/simil_matrix.determinant() -1.);
        // return rest_face_areas[f] * log((simil_matrix.transpose() * simil_matrix).trace()*0.5 / simil_matrix.determinant()); // 
        // return rest_face_areas[f] * current_triangle_area  * log((simil_matrix.transpose() * simil_matrix).trace()*0.5 / simil_matrix.determinant()); // 
        // return rest_face_areas[f] * current_triangle_area * (simil_matrix.transpose() * simil_matrix).trace()*0.5 / simil_matrix.determinant();
        
        
        // return rest_face_areas[f] * 
        //         (membrane_mu * simil_matrix.trace()/2. + membrane_lambda * simil_matrix.determinant()/4. -
        //          (membrane_mu/4. + membrane_lambda/2.) * log(simil_matrix.determinant()));
    });
    return membraneEnergy_func;
}


auto DeformationSolver::get_tinyAD_barrier_function(){
    auto barrierEnergy_func = TinyAD::scalar_function<3>(mesh->vertices()); // 
    barrierEnergy_func.add_elements<1>(mesh->vertices(), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {
        // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);

        // Get variables 3D vertex positions
        Vertex v = element.handle;
        Eigen::VectorX<T> p = element.variables(v);
        Eigen::VectorX<T> diffs = constraint_rhs - constraint_matrix * p;
        return -(diffs.array().log().sum());
    });
    return barrierEnergy_func;
}


void DeformationSolver::print_energies_after_transform(Eigen::Matrix3d A){
    // VertexData<Vector3> transformed_points(*mesh);
    DenseMatrix<double> mat_transformed_points = vertex_data_to_matrix(old_geometry->inputVertexPositions) * A.transpose();
    // for (Vertex v: mesh->vertices()){
    //     transformed_points[v] = vec_to_GC_vec3(mat_transformed_points.row(v.getIndex()));
    // }
    update_bending_rest_constants();
    update_membrane_rest_constants();
    // bending func tinyAD
    auto bendingEnergy_func = get_tinyAD_bending_function();
    // membrane func tinyAD
    auto membraneEnergy_func = get_tinyAD_membrane_function();

    Eigen::VectorXd x = tinyAD_flatten(mat_transformed_points);
    auto [membrane_f, membrane_g, membrane_H_proj] = membraneEnergy_func.eval_with_hessian_proj(x); //
    auto [bending_f, bending_g, bending_H_proj] = bendingEnergy_func.eval_with_hessian_proj(x); //
    
    TINYAD_DEBUG_OUT("\n\t\t membrane energy = " << membrane_f <<
                     ": bending= " << bending_f);
}


DenseMatrix<double> DeformationSolver::solve_for_bending(int visual_per_step, bool energy_plot, int* current_iter, float** energies_log){
    size_t num_var = 3 * mesh->nVertices();
    old_geometry->refreshQuantities();
    old_geometry->requireEdgeLengths();
    double initial_mean_edge_len = old_geometry->edgeLengths.toVector().mean();
    old_geometry->unrequireEdgeLengths();

    Eigen::SparseMatrix<double> cv_gauss_curve_diag;
    if (curvature_weighted_CP){
        convex_geometry->requireVertexGaussianCurvatures();
        Eigen::VectorXd convex_gaussian_curvatures = convex_geometry->vertexGaussianCurvatures.toVector();
        Eigen::VectorXd convex_gaussian_curvatures_for_flat(3*convex_gaussian_curvatures.size());
        convex_gaussian_curvatures_for_flat << convex_gaussian_curvatures, convex_gaussian_curvatures, convex_gaussian_curvatures; 
        std::cout << "convex_gaussian_curvatures min and max " << convex_gaussian_curvatures.minCoeff() << " " << convex_gaussian_curvatures.maxCoeff() << std::endl;
        cv_gauss_curve_diag = convex_gaussian_curvatures_for_flat.asDiagonal();
    }

    // tinyAD stuff
    build_constraint_matrix_and_rhs();
    auto barrierEnergy_func = get_tinyAD_barrier_function();
    update_bending_rest_constants();
    auto bendingEnergy_func = get_tinyAD_bending_function();
    update_membrane_rest_constants();
    auto membraneEnergy_func = get_tinyAD_membrane_function();    
    
    printf(" initializing variables\n");
    int n = mesh->nVertices();
    Eigen::VectorXd x = bendingEnergy_func.x_from_data([&] (Vertex v) {
        return to_eigen(old_geometry->inputVertexPositions[v.getIndex()]);
    });

    // GRBEnv env = GRBEnv();
    // GRBModel QPmodel = GRBModel(env);
    // build_QP_model_with_constraints(QPmodel, x, constraint_matrix, constraint_rhs);
            
    while (!check_feasibility(unflat_tinyAD(x))){ // assuming centered; which is
        printf("  -- scaling for feasibility -- \n");
        x *= 0.95;
    }

    VertexPositionGeometry *tmp_geometry = new VertexPositionGeometry(*mesh, unflat_tinyAD(x));
    
    // some parameters
    double barrier_lambda = init_barrier_lambda,
           bending_lambda = init_bending_lambda,
           membrane_lambda = init_membrane_lambda,
           CP_lambda = init_CP_lambda;
    internal_pt = 1.;

    vertex_only_assignment = VertexData<bool>(*convex_mesh, false);
    frozen_assignment = VertexData<Vertex>(*convex_mesh, Vertex());
    for (int i = 0; i < filling_max_iter; ++i) {
        // assign closest points
        std::vector<Vector3> post_frozen_points, pre_frozen_points;
        assign_closest_points_barycentric(tmp_geometry);
        for (Vertex cv: convex_mesh->vertices()){
            if (frozen_assignment[cv].getIndex() != INVALID_IND){
                pre_frozen_points.push_back(tmp_geometry->inputVertexPositions[frozen_assignment[cv]]);
                tmp_geometry->inputVertexPositions[frozen_assignment[cv]] = convex_geometry->inputVertexPositions[cv];
                post_frozen_points.push_back(convex_geometry->inputVertexPositions[cv]);
            }
        }
        if (enforce_snapping){
            polyscope::registerPointCloud("pre freeze points", pre_frozen_points)->setPointColor({1.,0.,0.})->setEnabled(true);
            polyscope::registerPointCloud("post freeze points", post_frozen_points)->setPointColor({0.,0.,1.})->setEnabled(true);
        }
        
        x = bendingEnergy_func.x_from_data([&] (Vertex v) {return to_eigen(tmp_geometry->inputVertexPositions[v.getIndex()]);});
        // polyscope::show();
        
        // DEBUG & visuals
        visualize_cv_cp_assignments(mesh, old_geometry, tmp_geometry, convex_mesh, convex_geometry, closest_point_assignment, closest_point_flat_operator, x);
        
        if (visual_per_step != 0){
            if (i % visual_per_step == 0){
                printf(" visualizing step %d, mesh size: %d\n", i, mesh->nVertices());
                auto tmp_PSmesh = polyscope::registerSurfaceMesh("temp sol", tmp_geometry->inputVertexPositions, mesh->getFaceVertexList());
                printf(" tmp geo container size %d\n", tmp_geometry->inputVertexPositions.size());
                printf(" ps mesh size %d\n", tmp_PSmesh->vertexDataSize);

                tmp_PSmesh->setSurfaceColor({136./255., 229./255., 107./255.});
                tmp_PSmesh->setEdgeWidth(1.);
                tmp_PSmesh->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
                tmp_PSmesh->setEnabled(true);
                // polyscope::screenshot();
                polyscope::frameTick();
            }
        }

        double CP_energy = closest_point_energy(tmp_geometry);
        // dynamic remesh 
        bool topo_changed = false;
        if (dynamic_remesh){
            // joint_remesh(mesh, old_geometry, tmp_geometry, initial_mean_edge_len);
            topo_changed = topo_changed || split_only_remesh(mesh, old_geometry, tmp_geometry, initial_mean_edge_len);
        }
        // split for small CP
        if (threshold_exceeded > 0){
            printf(" spliting\n");
            // modifies mesh, old_geometry, tmp_geometry
            topo_changed = topo_changed || split_barycentric_closest_points(tmp_geometry); 
        }
        // split when close to convex hull
        if (topo_changed){
            // update tineAD stuff
            barrierEnergy_func = get_tinyAD_barrier_function();
            update_bending_rest_constants();
            bendingEnergy_func = get_tinyAD_bending_function();
            update_membrane_rest_constants();
            membraneEnergy_func = get_tinyAD_membrane_function();
            x = bendingEnergy_func.x_from_data([&] (Vertex v) { // resized so need to reassign
                return to_eigen(tmp_geometry->inputVertexPositions[v.getIndex()]);
            });

            tmp_geometry = new VertexPositionGeometry(*mesh, unflat_tinyAD(x)); // TODO: for size reasons
            
            // re-assign CP since connectivity has changed
            assign_closest_points_barycentric(tmp_geometry);
            // assign_closest_vertices(tmp_geometry, true);
        }

        // CP stuff; splits shouldnt affect energy
        CP_energy = closest_point_energy(tmp_geometry);
        Eigen::VectorXd x_cv_flat = tinyAD_flatten(vertex_data_to_matrix(convex_geometry->inputVertexPositions));
        Eigen::SparseMatrix<double> A_CP = closest_point_flat_operator;
        if (curvature_weighted_CP){
            x_cv_flat = cv_gauss_curve_diag * x_cv_flat;
            A_CP = cv_gauss_curve_diag * A_CP;
        }
        Eigen::VectorX<double> CP_g = 2. * A_CP.transpose() * (A_CP * x - x_cv_flat);  // tinyAD_flatten(closest_point_energy_gradient(tmp_geometry));
        Eigen::SparseMatrix<double> CP_H = 2. * A_CP.transpose() * A_CP;


        // x - x_0 regularizer
        // Eigen::VectorXd reg_g = Eigen::VectorXd::Zero(x.size());
        Eigen::SparseMatrix<double> reg_H = 2.* identityMatrix<double>(x.size());
        

        // elastic stuff
        auto [bending_f, bending_g, bending_H_proj] = bendingEnergy_func.eval_with_hessian_proj(x); //
        auto [membrane_f, membrane_g, membrane_H_proj] = membraneEnergy_func.eval_with_hessian_proj(x); //
        // DEBUG?
        // auto [bending_f, bending_g, bending_H_proj] = bendingEnergy_func.eval_with_derivatives(x); //
        // auto [membrane_f, membrane_g, membrane_H_proj] = membraneEnergy_func.eval_with_derivatives(x); //
        
        // printf(" finding Barrier terms\n");
        // DenseMatrix<double> x_in_dense_format = vertex_data_to_matrix(tmp_geometry->inputVertexPositions);
        // auto [my_barrier_f, barrier_grad, barrier_hessian] = get_log_barrier_stuff(x_in_dense_format, Vector<double>::Ones(mesh->nVertices())); // thank u new c++
        // Eigen::SparseMatrix<double> my_barrier_H = tinyADify_barrier_hess(barrier_hessian);
        // Eigen::VectorX<double> my_barrier_g = tinyAD_flatten(barrier_grad);
          
        Eigen::SparseMatrix<double> total_H = bending_lambda * bending_H_proj + membrane_lambda * membrane_H_proj + CP_lambda * CP_H + reg_lambda * reg_H; // + barrier_lambda * barrier_H;
        Eigen::VectorXd total_g = bending_lambda * bending_g + membrane_lambda * membrane_g + CP_lambda * CP_g; // reg g is zero at x0 // + barrier_lambda * barrier_g;
        double total_energy = bending_lambda * bending_f + membrane_lambda * membrane_f + CP_lambda * CP_energy; // reg e is zero at x0 // barrier_lambda * barrier_f;
        // barrier stuff
        if (barrier_lambda != 0){
            auto [barrier_f, barrier_g, barrier_H] = barrierEnergy_func.eval_with_derivatives(x); //
            total_energy += barrier_lambda * barrier_f;
            total_g += barrier_lambda * barrier_g;
            total_H += barrier_lambda * barrier_H;
        }

        // TINYAD_INFO(" my barr vs TinyAD" << 
        //             "\n\t\t\tenergies: " << my_barrier_f << ", " << barrier_f <<
        //             "\n\t\t\tgrads: " << my_barrier_g.isApprox(barrier_g) << " "<< (barrier_g.transpose()-my_barrier_g.transpose()).norm()<<
        //             "\n\t\t\thessians:" << my_barrier_H.isApprox(barrier_H) << " "<< (barrier_H - my_barrier_H).norm());

        TINYAD_DEBUG_OUT("\t- Energy in iter " << i << ": bending= " << bending_lambda << "*" << bending_f << 
                         "\n\t\t\t\t membrane  = " << membrane_lambda << " * " << membrane_f <<
                         "\n\t\t\t\t CP        = " << CP_lambda << " * "<< CP_energy);
        if (barrier_lambda != 0.) TINYAD_DEBUG_OUT("\n\t\t\t\t barrier   = " << barrier_lambda << " * " << barrierEnergy_func.eval(x));
        TINYAD_DEBUG_OUT("\n\t\t\t\t\t total: " << bending_lambda * bending_f + membrane_lambda * membrane_f + CP_lambda * CP_energy);
        if (energy_plot){
            *current_iter = i;
            energies_log[0][i] = bending_f;
            energies_log[1][i] = membrane_f;
            energies_log[2][i] = CP_energy;
        }
        Eigen::VectorXd old_x = x;
        Eigen::VectorXd d;
        if (barrier_lambda == 0){ // use QP
            // Eigen::VectorXd new_x = update_QP_objective_and_solve(QPmodel, total_H/2., total_g - total_H * old_x, old_x);
            Eigen::MatrixX<bool> active_set = get_active_set_matrix(unflat_tinyAD(x), active_set_threshold);
            // std::cout << ANSI_FG_GREEN << "total active consts: " << active_set.cast<int>().sum() << "/" << n*constraint_matrix.rows() << ANSI_RESET << std::endl; 
            Eigen::VectorXd new_x = solve_QP_with_ineq_GRB(total_H/2., total_g - total_H * x, old_x, //
                                                       constraint_matrix, constraint_rhs, 
                                                       get_frozen_flags(), get_frozen_x(),
                                                       active_set);
            d = new_x - old_x;
        }
        else { // my own barrier
            PositiveDefiniteSolver<double> newtonSolver(total_H);
            d = newtonSolver.solve(-total_g);
        }
        x = line_search(old_x, d, total_energy, total_g,
                        [&] (const Eigen::VectorXd curr_x, bool print = false) {
                                if (!check_feasibility(unflat_tinyAD(curr_x)))
                                    return std::numeric_limits<double>::infinity();
                                double bending_e = bendingEnergy_func.eval(curr_x),
                                    membrane_e = membraneEnergy_func.eval(curr_x),
                                    CP_e = closest_point_energy(curr_x),
                                    reg_e = (curr_x - old_x).squaredNorm();
                                double f_new = bending_lambda * bending_e + membrane_lambda * membrane_e + CP_lambda * CP_e + reg_lambda * reg_e; // + barrier_lambda * log_barr_e;  
                                if (barrier_lambda != 0.)
                                    f_new += barrier_lambda * barrierEnergy_func.eval(curr_x);
                                return f_new;
                            },
                        1., 0.9, 200, 0.);
                    // smax, decay, max_iter, armijo constant
        // update tmp geometry
        bendingEnergy_func.x_to_data(x, [&] (Vertex v, const Eigen::Vector3d& p) {
            tmp_geometry->inputVertexPositions[v] = to_geometrycentral(p);
        });

        double step_norm = (x - old_x).norm();
        // printf(" step norm is %9f\n", step_norm);
        internal_pt *= internal_growth_p;
        bending_lambda = get_scheduled_weight(init_bending_lambda, final_bending_lambda);
        membrane_lambda = get_scheduled_weight(init_membrane_lambda, final_membrane_lambda);
        CP_lambda = get_scheduled_weight(init_CP_lambda, final_CP_lambda);
        barrier_lambda = get_scheduled_weight(init_barrier_lambda, final_barrier_lambda);
        // if (step_norm < convergence_eps)
        //     break;
        if (step_norm == 0 && bending_lambda == final_bending_lambda && membrane_lambda == final_membrane_lambda && CP_lambda == final_CP_lambda && barrier_lambda == final_barrier_lambda)
            break;
    }
    // last tick
    auto tmp_PSmesh = polyscope::registerSurfaceMesh("temp sol", tmp_geometry->inputVertexPositions, mesh->getFaceVertexList());
    tmp_PSmesh->setSurfaceColor({136./255., 229./255., 107./255.});
    tmp_PSmesh->setEdgeWidth(1.);
    tmp_PSmesh->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
    tmp_PSmesh->setEnabled(true);
    // polyscope::screenshot();
    polyscope::frameTick();

    DenseMatrix<double> new_points_mat = unflat_tinyAD(x); 
    deformed_geometry = new VertexPositionGeometry(*mesh, new_points_mat); // TODO: update earlier? per temp geo?
    return new_points_mat;
}



VertexData<DenseMatrix<double>> DeformationSolver::per_vertex_G_jacobian(VertexPositionGeometry *tmp_geometry){
    DenseMatrix<double> zmat = DenseMatrix<double>::Zero(3,3);
    VertexData<DenseMatrix<double>> dG_dv(*mesh, zmat);
    for (Face f: mesh->faces()){
        double face_area = tmp_geometry->faceArea(f);
        Vector3 face_normal = tmp_geometry->faceNormal(f); // assuming outward normals

        size_t face_degree = f.degree();
        Vector3 vert_sum = Vector3::zero(); 
        for (Vertex tmp_v: f.adjacentVertices())
            vert_sum += tmp_geometry->inputVertexPositions[tmp_v];
        for (Halfedge he: f.adjacentHalfedges()){
            Vertex v = he.tailVertex();
            Vector3 p = tmp_geometry->inputVertexPositions[v];
            Vector3 Gf_G = (vert_sum + p)/(double)(face_degree + 1) - current_G;
            DenseMatrix<double> tmp_mat = vec32vec(Gf_G) * 
                                          vec32vec(face_normal).transpose();
            assert(tmp_mat.cols() == 3);
            assert(tmp_mat.rows() == 3);
            dG_dv[v] += face_area * tmp_mat;
        }
    }
    for (Vertex v: mesh->vertices())
        dG_dv[v] /= 3.*current_volume;
    return dG_dv;
}

Eigen::VectorXd DeformationSolver::flat_distance_multiplier(VertexPositionGeometry *tmp_geometry, bool from_faces){
    size_t n = mesh->nVertices();
    Eigen::VectorXd flat_dist_vec = Eigen::VectorXd::Zero(3*n);
    VertexData<double> dists(*mesh, 0.);
    for (Vertex v: mesh->vertices()){
        double distance_to_cv = get_point_distance_to_convex_hull(vec32vec(tmp_geometry->inputVertexPositions[v]), from_faces);
        dists[v] = distance_to_cv;
        flat_dist_vec[3 * v.getIndex() + 0] = distance_to_cv;
        flat_dist_vec[3 * v.getIndex() + 1] = distance_to_cv;
        flat_dist_vec[3 * v.getIndex() + 2] = distance_to_cv;
    }
    // polyscope::getSurfaceMesh("temp sol")->addVertexScalarQuantity("dist to CV", dists)->setEnabled(true);
    return flat_dist_vec;
}

Eigen::MatrixXd DeformationSolver::per_vertex_G_derivative(VertexPositionGeometry *tmp_geometry){
    size_t n = mesh->nVertices();

    VertexData<DenseMatrix<double>> dG_dv = per_vertex_G_jacobian(tmp_geometry);
    Eigen::MatrixXd G_Ghat_dV = Eigen::MatrixXd::Zero(n,3);

    Vector3 goal_G_dir = current_G - goal_G;
    for (Vertex v: mesh->vertices()){
        G_Ghat_dV.row(v.getIndex()) = 2. * dG_dv[v].transpose() * vec32vec(goal_G_dir);
    }
    return G_Ghat_dV;
}


DenseMatrix<double> DeformationSolver::solve_for_G(int visual_per_step, 
                                                   bool energy_plot, int* current_iter, float** energies_log){
    std::cout << "solving for G: \n -- mesh size "<< mesh->nVertices() << "\n";
    
    // some parameters
    size_t num_var = 3 * mesh->nVertices();
    old_geometry->refreshQuantities();
    old_geometry->requireEdgeLengths();
    double initial_mean_edge_len = old_geometry->edgeLengths.toVector().mean();
    old_geometry->unrequireEdgeLengths();

    // tinyAD stuff
    build_constraint_matrix_and_rhs();
    update_bending_rest_constants();
    auto bendingEnergy_func = get_tinyAD_bending_function();
    update_membrane_rest_constants();
    auto membraneEnergy_func = get_tinyAD_membrane_function();    
    
    printf(" initializing variables\n");
    int n = mesh->nVertices();
    Eigen::VectorXd x = bendingEnergy_func.x_from_data([&] (Vertex v) {
        return to_eigen(deformed_geometry->inputVertexPositions[v.getIndex()]);
    });
    VertexPositionGeometry *tmp_geometry = new VertexPositionGeometry(*mesh, unflat_tinyAD(x));
    assert(check_feasibility(unflat_tinyAD(x))); // the deformed shape should be feasible
    std::cout << "temp geo vs def geo: \n -- "<< tmp_geometry->mesh.nVertices() << "\n -- " << deformed_geometry->mesh.nVertices() << "\n";

    // some parameters
    double bending_lambda = final_bending_lambda,
           membrane_lambda = final_membrane_lambda,
           G_lambda = init_G_lambda;
    internal_pt = 1.;

    for (int i = 0; i < filling_max_iter; ++i) {    
        // update G
        auto GV_pair = find_center_of_mass(*mesh, *tmp_geometry);
        current_G = GV_pair.first;
        current_volume = GV_pair.second;
        // visuals
        if (visual_per_step != 0){
            if ((i+1) % visual_per_step == 0){
                printf(" visualizing step %d\n", i);
                polyscope::registerPointCloud("G", std::vector<Vector3>({current_G}))->setPointColor({1.,0.,0.})->setPointRadius(0.01)->setEnabled(true);
                polyscope::registerPointCloud("goal G", std::vector<Vector3>({goal_G}))->setPointColor({0.,0.,1.})->setPointRadius(0.01)->setEnabled(true);
                
                auto tmp_PSmesh = polyscope::registerSurfaceMesh("temp sol", tmp_geometry->inputVertexPositions, mesh->getFaceVertexList());
                tmp_PSmesh->setSurfaceColor({136./255., 229./255., 107./255.});
                tmp_PSmesh->setEdgeWidth(1.);
                tmp_PSmesh->setTransparency(0.7);
                tmp_PSmesh->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
                tmp_PSmesh->setEnabled(true);
                // polyscope::screenshot();
                polyscope::frameTick();
            }
        }

        bool topo_changed = false;
        // if (dynamic_remesh){
        //     // joint_remesh(mesh, old_geometry, tmp_geometry, initial_mean_edge_len);
        //     topo_changed = topo_changed || split_only_remesh(mesh, old_geometry, tmp_geometry, initial_mean_edge_len);
        // }
        // if (topo_changed){
        //     std::cout << ANSI_FG_RED << "topo changed" << ANSI_RESET << std::endl;
        //     // update tineAD stuff
        //     update_bending_rest_constants();
        //     bendingEnergy_func = get_tinyAD_bending_function();
        //     update_membrane_rest_constants();
        //     membraneEnergy_func = get_tinyAD_membrane_function();
        //     x = bendingEnergy_func.x_from_data([&] (Vertex v) { // resized so need to reassign
        //         return to_eigen(tmp_geometry->inputVertexPositions[v.getIndex()]);
        //     });
        // }
        // get G diff energy and derivative
        double Gdiff_f = (current_G - goal_G).norm2();
        Eigen::VectorXd Gdiff_g = tinyAD_flatten(per_vertex_G_derivative(tmp_geometry));

        // ||x - x_0|| regularizer
        // Eigen::VectorXd reg_g = Eigen::VectorXd::Zero(x.size());
        Eigen::SparseMatrix<double> reg_H = 2.* identityMatrix<double>(x.size());
        
        // elastic stuff
        auto [bending_f, bending_g] = bendingEnergy_func.eval_with_gradient(x); //
        auto [membrane_f, membrane_g] = membraneEnergy_func.eval_with_gradient(x); //

        Eigen::VectorXd total_g = bending_lambda * bending_g + membrane_lambda * membrane_g + G_lambda * Gdiff_g; // reg g is zero at x0 // + barrier_lambda * barrier_g;
        double total_energy = bending_lambda * bending_f + membrane_lambda * membrane_f + G_lambda * Gdiff_f; // reg e is zero at x0 // barrier_lambda * barrier_f;

        Eigen::VectorXd flat_dist_mult = flat_distance_multiplier(tmp_geometry, true);
        
        flat_dist_mult /= flat_dist_mult.maxCoeff();
        assert(flat_dist_mult.size() == x.size());
        // DEBUG
        // printf("temp sol size %d\n other shit size", polyscope::getSurfaceMesh("temp sol")->vertexDataSize, unflat_tinyAD(flat_dist_mult).col(0).size());
        polyscope::getSurfaceMesh("temp sol")->addVertexScalarQuantity("normalized distance to CV", unflat_tinyAD(flat_dist_mult).col(0))->setEnabled(false);
        polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("total g un-wieghted", unflat_tinyAD(total_g))->setEnabled(false);
        
        total_g = total_g.cwiseProduct(flat_dist_mult.cwiseAbs2());

        // std::cout << "temp sol size:" << mesh->nVertices() << " gghat size: " << G_Ghat_g.size() << "unflat gghat size" << unflat_tinyAD(G_Ghat_g).rows() << std::endl;
        polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("bending g", unflat_tinyAD(bending_g))->setEnabled(false);
        polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("membrane g", unflat_tinyAD(membrane_g))->setEnabled(false);
        polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("Gdiff un-wieghted", unflat_tinyAD(Gdiff_g))->setEnabled(true);
        polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("total g wieghted", unflat_tinyAD(total_g))->setEnabled(true);
        polyscope::frameTick();

        TINYAD_DEBUG_OUT("\t- Energy in iter " << i << ": bending= " << bending_lambda << "*" << bending_f << 
                         "\n\t\t\t\t membrane  = " << membrane_lambda << " * " << membrane_f <<
                         "\n\t\t\t\t G diff    = " << G_lambda << " * "<< Gdiff_f);
        TINYAD_DEBUG_OUT("\n\t\t\t\t\t total: " << bending_lambda * bending_f + membrane_lambda * membrane_f + G_lambda * Gdiff_f);
        if (energy_plot){
            *current_iter = i;
            energies_log[0][i] = bending_f;
            energies_log[1][i] = membrane_f;
            energies_log[2][i] = Gdiff_f;
        }
        Eigen::VectorXd old_x = x;
        Eigen::VectorXd d = -total_g;
        x = line_search(old_x, d, total_energy, total_g,
                        [&] (const Eigen::VectorXd curr_x, bool print = false) {
                                if (!check_feasibility(unflat_tinyAD(curr_x)))
                                    return std::numeric_limits<double>::infinity();
                                // find new G at every step
                                VertexPositionGeometry *tmp_tmp_geo = new VertexPositionGeometry(*mesh, unflat_tinyAD(curr_x));
                                Vector3 new_G = find_center_of_mass(*mesh, *tmp_tmp_geo).first;
                                
                                double bending_e = bendingEnergy_func.eval(curr_x),
                                       membrane_e = membraneEnergy_func.eval(curr_x),
                                       G_Ghat_e = (new_G - goal_G).norm2(),
                                       reg_e = (curr_x - old_x).squaredNorm();
                                double f_new = bending_lambda * bending_e + membrane_lambda * membrane_e + G_lambda * G_Ghat_e + reg_lambda * reg_e; // + barrier_lambda * log_barr_e;  
                                return f_new;
                            },
                        1., 0.9, 200, 0.);
                    // smax, decay, max_iter, armijo constant
        // update tmp geometry
        bendingEnergy_func.x_to_data(x, [&] (Vertex v, const Eigen::Vector3d& p) {
            tmp_geometry->inputVertexPositions[v] = to_geometrycentral(p);
        });

        internal_pt *= internal_growth_p;
        G_lambda = get_scheduled_weight(init_G_lambda, final_G_lambda);
        
        double step_norm = (x - old_x).norm();
        // if (step_norm == 0. && G_lambda == final_G_lambda)
        //     break;
    }
    DenseMatrix<double> new_points_mat = unflat_tinyAD(x); 
    return new_points_mat;
}




/// Visuals


void visualize_cv_cp_assignments(ManifoldSurfaceMesh *mesh, 
                                 VertexPositionGeometry *old_geometry, VertexPositionGeometry *tmp_geometry, 
                                 ManifoldSurfaceMesh *convex_mesh, VertexPositionGeometry *convex_geometry, 
                                 VertexData<SurfacePoint> closest_point_assignment, 
                                 Eigen::SparseMatrix<double> closest_point_flat_operator, Eigen::VectorXd x){
    
    polyscope::registerPointCloud("CVs", convex_geometry->inputVertexPositions)->setPointColor({1.,0.,0.})->setEnabled(false);
    polyscope::registerPointCloud("CPs", unflat_tinyAD(closest_point_flat_operator * x))->setPointColor({0.,1.,0.})->setEnabled(false);
    
    std::vector<std::array<size_t, 2>> edge_inds;
    std::vector<Vector3> edge_points;
    std::vector<double> edge_lens;
    for (Vertex v: convex_mesh->vertices()){
        edge_inds.push_back({2*v.getIndex(),2*v.getIndex() + 1});
        edge_points.push_back(convex_geometry->inputVertexPositions[v]);
        edge_points.push_back(closest_point_assignment[v].interpolate(tmp_geometry->inputVertexPositions));
        edge_lens.push_back((edge_points.back() - edge_points[edge_points.size()-2]).norm());
    }
    polyscope::registerCurveNetwork("CV-CP", edge_points, edge_inds)->setColor({0.,0.,1.})->setRadius(0.0003)->addEdgeScalarQuantity("edge len", edge_lens);
    polyscope::show();
}