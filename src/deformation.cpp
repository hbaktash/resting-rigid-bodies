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

#include "deformation.h"

#include "geometrycentral/numerical/linear_solvers.h"

// Gradient stuff
// formula source in header file
Vector3 dihedral_angle_grad_G(Vector3 G, Vector3 A, Vector3 B, Vector3 C){
    // A is in the middle (on the unit sphere)
    Vector3 GA = A - G,
            GB = B - G,
            GC = C - G;
    // un-normalized normals
    // order: B A C
    Vector3 nB_un = cross(GA, B - A),
            nC_un = cross(GA, A - C);
    Vector3 grad1 = nB_un * dot(-GA, B - A)/nB_un.norm2(),
            grad2 = nC_un * dot(GA, A - C)/nC_un.norm2();
    return -(grad1 + grad2)/GA.norm();
}
// same as G, replacing G and A
Vector3 dihedral_angle_grad_A(Vector3 G, Vector3 A, Vector3 B, Vector3 C){
    return -dihedral_angle_grad_G(A, G, B, C);
}
// wing side vertex B; kinda simple
Vector3 dihedral_angle_grad_B(Vector3 G, Vector3 A, Vector3 B, Vector3 C){
    Vector3 nB_un = cross(A - G, B - A);
    return  nB_un * (G-A).norm()/nB_un.norm2();
}
// wing side vertex C; kinda simple
Vector3 dihedral_angle_grad_C(Vector3 G, Vector3 A, Vector3 B, Vector3 C){
    Vector3 nC_un = cross(A - G, A - C);
    return nC_un * (G-A).norm()/nC_un.norm2();
}


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

geometrycentral::DenseMatrix<double> get_ARAP_positions(
                                       geometrycentral::DenseMatrix<double> old_pos_mat,
                                       geometrycentral::DenseMatrix<double> new_pos_mat, 
                                       geometrycentral::DenseMatrix<double> init_sol, 
                                       ManifoldSurfaceMesh &inner_mesh,
                                       geometrycentral::Vector<int> hull_indices){
    Eigen::MatrixXd V,U, bc;
    Eigen::MatrixXi F;
    Eigen::VectorXi S,b;
    igl::ARAPData arap_data;
    arap_data.max_iter = 20;
    SparseMatrix<double> L;
    F = inner_mesh.getFaceVertexMatrix<size_t>().cast<int>();
    V = old_pos_mat;
    b = hull_indices;
    // igl::arap_precomputation(old_pos_mat, F, 3, hull_indices, arap_data);
    igl::arap_precomputation(V, F, 3, b, arap_data);
    // geometrycentral::DenseMatrix<double> U = init_sol; 
    U = init_sol;
    bc = new_pos_mat(hull_indices, Eigen::all);
    igl::arap_solve(bc, arap_data, U);
    return U;
}


DeformationSolver::DeformationSolver(ManifoldSurfaceMesh *_mesh, VertexPositionGeometry *_old_geometry,
                                     ManifoldSurfaceMesh *_convex_mesh, VertexPositionGeometry *_convex_geometry){
    mesh = _mesh;
    old_geometry = _old_geometry;
    convex_mesh = _convex_mesh;
    convex_geometry = _convex_geometry;
}

double DeformationSolver::get_scheduled_weight(double init_w, double final_w){
    return final_w + (init_w - final_w) * internal_pt * internal_growth_p;
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
    printf(" CP energy: %f \n", energy);
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
        
        // update assignment and involvement
        closest_point_assignment[c_v] = SurfacePoint(closest_vertex); 
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
    closest_point_assignment = VertexData<SurfacePoint>(*convex_mesh, SurfacePoint());
    CP_involvement = VertexData<bool>(*mesh, false);
    // stats
    face_assignments = 0; edge_assignments = 0; vertex_assignments = 0;

    // linear operator corresponding to assignments
    closest_point_operator = Eigen::SparseMatrix<double>(convex_mesh->nVertices(), mesh->nVertices());
    closest_point_flat_operator = Eigen::SparseMatrix<double>(3 * convex_mesh->nVertices(), 3*mesh->nVertices());
    std::vector<Eigen::Triplet<double>> tripletList, flat_tripletList;
    tripletList.reserve(3*mesh->nVertices());
    flat_tripletList.reserve(9*mesh->nVertices());

    for (Vertex c_v: convex_mesh->vertices()){
        Face closest_face;
        Edge closest_edge;
        Vertex closest_vertex;
        
        Vector3 on_edge_projection;
        double min_face_dist  = std::numeric_limits<double>::infinity(),
               min_edge_dist   = std::numeric_limits<double>::infinity(),
               min_vertex_dist  = std::numeric_limits<double>::infinity();

        Vector3 c_p = convex_geometry->inputVertexPositions[c_v];

        // there should be a smarter way of checking all elements in one loop; maybe even without a flagging them??
        for (Face f: mesh->faces()){ 
            Vector3 f_normal = new_geometry->faceNormal(f);
            Halfedge curr_he  = f.halfedge(),
                    first_he = f.halfedge();
            // assume outward normals; which is inward w.r.t. p
            double face_dist = dot(f_normal, c_p - new_geometry->inputVertexPositions[f.halfedge().vertex()]); // using some point of f
            
            bool face_is_projectable = true;
            while (true){ // checking if face-projectable
                Vertex v1 = curr_he.tailVertex(), v2 = curr_he.tipVertex();
                Vector3 A = new_geometry->inputVertexPositions[v1], 
                        B = new_geometry->inputVertexPositions[v2];
                Vector3 N_PAB = cross(B - A, A - c_p);
                if (dot(f_normal, N_PAB) <= 0) { // not face-projectable on this face
                    face_is_projectable = false;
                    break; // go to next face
                }
                curr_he = curr_he.next();
                if (curr_he == first_he)
                    break;
            }
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
        if (min_dist == min_face_dist){ // assuming triangle mesh
            // printf("closest is a face\n");
            face_assignments++;

            Vector3 f_normal = new_geometry->faceNormal(closest_face);
            Vector3 on_face_projection = c_p - f_normal * dot(f_normal, 
                                                              c_p - new_geometry->inputVertexPositions[closest_face.halfedge().vertex()]);
            Vertex v1 = closest_face.halfedge().vertex(),
                   v2 = closest_face.halfedge().next().vertex(),
                   v3 = closest_face.halfedge().next().next().vertex();
            Vector3 A = new_geometry->inputVertexPositions[v1],
                    B = new_geometry->inputVertexPositions[v2],
                    C = new_geometry->inputVertexPositions[v3];
            // SurfacePoint assignment
            Vector3 bary_coor = barycentric(on_face_projection, A, B, C);
            closest_point_assignment[c_v] = SurfacePoint(closest_face, bary_coor);

            // update CP involvement
            CP_involvement[v1] = true;
            CP_involvement[v2] = true;
            CP_involvement[v3] = true;

            // operator entries
            tripletList.emplace_back(c_v.getIndex(), v1.getIndex(), bary_coor.x);
            tripletList.emplace_back(c_v.getIndex(), v2.getIndex(), bary_coor.y);
            tripletList.emplace_back(c_v.getIndex(), v3.getIndex(), bary_coor.z);
            for (size_t d = 0; d < 3; d++){
                flat_tripletList.emplace_back(3 * c_v.getIndex() + d, 3*v1.getIndex() + d, bary_coor.x);
                flat_tripletList.emplace_back(3 * c_v.getIndex() + d, 3*v1.getIndex() + d, bary_coor.y);
                flat_tripletList.emplace_back(3 * c_v.getIndex() + d, 3*v1.getIndex() + d, bary_coor.z);
            }
        }
        else if (min_dist == min_edge_dist){
            // printf("closest is an edge\n");
            edge_assignments++;

            Vertex v1 = closest_edge.firstVertex(),
                   v2 = closest_edge.secondVertex();
            Vector3 A = new_geometry->inputVertexPositions[v1], 
                    B = new_geometry->inputVertexPositions[v2];
            double tVal = (on_edge_projection - A).norm()/(A - B).norm();
            closest_point_assignment[c_v] = SurfacePoint(closest_edge, tVal);

            // update CP involvement
            CP_involvement[v1] = true;
            CP_involvement[v2] = true;
            
            // operator entries
            tripletList.emplace_back(c_v.getIndex(), v1.getIndex(), 1. - tVal);
            tripletList.emplace_back(c_v.getIndex(), v2.getIndex(), tVal);
            for (size_t d = 0; d < 3; d++){
                flat_tripletList.emplace_back(3 * c_v.getIndex() + d, 3*v1.getIndex() + d, 1. - tVal);
                flat_tripletList.emplace_back(3 * c_v.getIndex() + d, 3*v2.getIndex() + d, tVal);
            }
        }
        else {// min_dist == min_vertex_dist
            // printf("closest is a vertex\n");
            vertex_assignments++;

            closest_point_assignment[c_v] = SurfacePoint(closest_vertex); 
            // update CP involvement
            CP_involvement[closest_vertex] = true;
              
            tripletList.emplace_back(c_v.getIndex(), closest_vertex.getIndex(), 1.);
            for (size_t d = 0; d < 3; d++)
                flat_tripletList.emplace_back(3 * c_v.getIndex() + d, 3*closest_vertex.getIndex() + d, 1.);
        }
    }
    closest_point_operator.setFromTriplets(tripletList.begin(), tripletList.end());
    closest_point_flat_operator.setFromTriplets(flat_tripletList.begin(), flat_tripletList.end());
    printf("  ** stats for CP: vertices %d, edges: %d, faces: %d \n", vertex_assignments, edge_assignments, face_assignments);
}


void DeformationSolver::split_barycentric_closest_points(VertexPositionGeometry *new_geometry){  
    // assing CP should be called already
    size_t mutli_edge_hits = 0, multi_face_hits = 0, multi_vertex_hits = 0;
    EdgeData<bool> edge_is_hit(*mesh, false);
    FaceData<bool> face_is_hit(*mesh, false);
    VertexData<bool> vertex_is_hit(*mesh, false);

    // auto CPpcold = polyscope::registerPointCloud("closest points on old geo", closest_point_operator*vertex_data_to_matrix(old_geometry->inputVertexPositions));
    // auto CPpcnew = polyscope::registerPointCloud("closest points on new geo", closest_point_operator*vertex_data_to_matrix(new_geometry->inputVertexPositions));
    // CPpcold->setPointRadius( 0.02, false);
    // CPpcold->setPointColor( {0.2,1,0});
    // CPpcnew->setPointRadius( 0.02, false);
    // CPpcnew->setPointColor( {0.2,0.1,1});
    
    for (Vertex cv: convex_mesh->vertices()){
        SurfacePoint assigned_cp = closest_point_assignment[cv];
        Vertex new_v;
        Vector3 split_p_old_geo,
                split_p_new_geo;
        if (assigned_cp.type == SurfacePointType::Edge){
            if (!edge_is_hit[assigned_cp.edge])
                edge_is_hit[assigned_cp.edge] = true;
            else {
                mutli_edge_hits++;
                continue;
            }
            // printf("splitting edge %d: %d,%d \n", assigned_cp.edge.getIndex(), assigned_cp.edge.firstVertex().getIndex(), assigned_cp.edge.secondVertex().getIndex());
            // must be done before spliting
            split_p_old_geo = assigned_cp.interpolate(old_geometry->inputVertexPositions);
            split_p_new_geo = assigned_cp.interpolate(new_geometry->inputVertexPositions);
            // debug 
            new_v = mesh->splitEdgeTriangular(assigned_cp.edge).vertex();

        }
        else if (assigned_cp.type == SurfacePointType::Face){
            if (!face_is_hit[assigned_cp.face])
                face_is_hit[assigned_cp.face] = true;
            else {
                multi_face_hits++;
                continue;
            } 
            // printf("splitting face %d\n", assigned_cp.face.getIndex());
            // must be done before spliting
            split_p_old_geo = assigned_cp.interpolate(old_geometry->inputVertexPositions);
            split_p_new_geo = assigned_cp.interpolate(new_geometry->inputVertexPositions);
            
            new_v = mesh->insertVertex(assigned_cp.face);
        }
        else {
            // printf("splitting vertex %d\n", assigned_cp.vertex.getIndex());
            // must be done before spliting
            split_p_old_geo = assigned_cp.interpolate(old_geometry->inputVertexPositions);
            split_p_new_geo = assigned_cp.interpolate(new_geometry->inputVertexPositions);
            new_v = assigned_cp.vertex;
            if (!vertex_is_hit[assigned_cp.vertex])
                vertex_is_hit[assigned_cp.vertex] = true;
            else 
                multi_vertex_hits++;
        }
        old_geometry->inputVertexPositions[new_v] = split_p_old_geo;
        new_geometry->inputVertexPositions[new_v] = split_p_new_geo;

        // std::cout<< "\t\t new p : pmid " << old_geometry->inputVertexPositions[new_v] <<"\n";
    
    }

    mesh->compress();
    // polyscope::registerSurfaceMesh("old geo post split", old_geometry->inputVertexPositions, mesh->getFaceVertexList());
    // polyscope::registerSurfaceMesh("new geo post split", new_geometry->inputVertexPositions, mesh->getFaceVertexList());
    // polyscope::show();
    printf(" splitting is over. multi edge hits %d, \n\t\t\t\t multi face hits: %d, \n\t\t\t\t multi vertex hits: %d\n", mutli_edge_hits, multi_face_hits, multi_vertex_hits);
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
            printf(" edge is dead %d\n", e.isDead());
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
        
        return rest_face_areas[f] * (simil_matrix.transpose() * simil_matrix).trace() / simil_matrix.determinant(); // 
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


DenseMatrix<double> DeformationSolver::solve_for_bending(int visual_per_step){
    size_t num_var = 3 * mesh->nVertices();
    old_geometry->refreshQuantities();
    old_geometry->requireEdgeLengths();
    double initial_mean_edge_len = old_geometry->edgeLengths.toVector().mean();
    old_geometry->unrequireEdgeLengths();

    // barrier func tinyAD
    build_constraint_matrix_and_rhs();
    auto barrierEnergy_func = get_tinyAD_barrier_function();

    // // bending rest constants
    update_bending_rest_constants();
    // bending func tinyAD
    auto bendingEnergy_func = get_tinyAD_bending_function();

    // // membrane rest constants
    update_membrane_rest_constants();
    // membrane func tinyAD
    auto membraneEnergy_func = get_tinyAD_membrane_function();
    
    
    printf(" initializing variables\n");
    int n = mesh->nVertices();
    Eigen::VectorXd x = bendingEnergy_func.x_from_data([&] (Vertex v) {
        return to_eigen(old_geometry->inputVertexPositions[v.getIndex()]);
    });

    while (!check_feasibility(unflat_tinyAD(x))){ // assuming centered; which is
        printf("  -- scaling for feasibility -- \n");
        x *= 0.95;
    }

    VertexPositionGeometry *tmp_geometry = new VertexPositionGeometry(*mesh);
    bendingEnergy_func.x_to_data(x, [&] (Vertex v, const Eigen::Vector3d& p) {
        tmp_geometry->inputVertexPositions[v] = to_geometrycentral(p);
    });

    
    // some parameters
    double convergence_eps = 1e-7;
    double barrier_lambda = init_barrier_lambda,
           bending_lambda = init_bending_lambda,
           membrane_lambda = init_membrane_lambda,
           CP_lambda = init_CP_lambda;
    bool do_barycentric = true;

    for (int i = 0; i < filling_max_iter; ++i) {
        
        // visuals
        if (visual_per_step != 0){
            if (i % visual_per_step == 0){
                // printf(" visualizing step %d\n", i);
                
                auto tmp_PSmesh = polyscope::registerSurfaceMesh("temp sol", tmp_geometry->inputVertexPositions, mesh->getFaceVertexList());
                tmp_PSmesh->setSurfaceColor({136./255., 229./255., 107./255.});
                tmp_PSmesh->setEdgeWidth(1.);
                tmp_PSmesh->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
                tmp_PSmesh->setEnabled(true);
                polyscope::frameTick();
            }
        }
        // assign closest points
        printf(" assigning  closest points\n");
        if (do_barycentric)
            assign_closest_points_barycentric(tmp_geometry);
        else 
            assign_closest_vertices(tmp_geometry);
        // printf(" CP energy\n");
        double CP_energy = closest_point_energy(tmp_geometry);

        // split when close to convex hull
        if (CP_energy <= refinement_CP_threshold && do_barycentric){
            polyscope::warning("low CP energy " + std::to_string(CP_energy), "Spliting SurfacePoints and switching to Vertex only assignment");
            printf(" spliting\n");
            // modifies mesh, old_geometry, tmp_geometry
            split_barycentric_closest_points(tmp_geometry); 
            do_barycentric = false;

            // barrier func tinyAD
            // convex frame has not changed
            barrierEnergy_func = get_tinyAD_barrier_function();

            // // bending rest constants
            update_bending_rest_constants();
            // bending func tinyAD
            bendingEnergy_func = get_tinyAD_bending_function();

            // // membrane rest constants
            update_membrane_rest_constants();
            // membrane func tinyAD
            membraneEnergy_func = get_tinyAD_membrane_function();

            // re-build x; since the size has changed 
            x = bendingEnergy_func.x_from_data([&] (Vertex v) {
                return to_eigen(tmp_geometry->inputVertexPositions[v.getIndex()]);
            });

            // re assign CP since connectivity has changed
            // assign_closest_points_barycentric(tmp_geometry);
            assign_closest_vertices(tmp_geometry, true);
        }
        
        // printf(" finding CP terms\n");

        // CP stuff; splits shouldnt affect energy
        CP_energy = closest_point_energy(tmp_geometry); // making sure subdivision has gone right
        Eigen::VectorXd x_cv_flat = tinyAD_flatten(vertex_data_to_matrix(convex_geometry->inputVertexPositions));
        Eigen::SparseMatrix<double> A_CP = closest_point_flat_operator;
        Eigen::VectorX<double> CP_g = 2. * A_CP.transpose() * (A_CP * x - x_cv_flat);  // tinyAD_flatten(closest_point_energy_gradient(tmp_geometry));
        Eigen::SparseMatrix<double> CP_H = 2. * A_CP.transpose() * A_CP;
        
        // x - x_0 regularizer
        // Eigen::VectorXd reg_g = Eigen::VectorXd::Zero(x.size());
        Eigen::SparseMatrix<double> reg_H = 2.* identityMatrix<double>(x.size());
        

        // polyscope::registerPointCloud("closest assignment from flat mat", unflat_tinyAD(A_CP * x));
        // polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("CP grads", -unflat_tinyAD(CP_g))->setEnabled(true);
        // polyscope::show();
        // bending stuff
        printf(" bending Hessian eval\n ");
        auto [bending_f, bending_g, bending_H_proj] = bendingEnergy_func.eval_with_hessian_proj(x); //
        printf(" membrane Hessian eval\n ");
        auto [membrane_f, membrane_g, membrane_H_proj] = membraneEnergy_func.eval_with_hessian_proj(x); //
        
        // barrier stuff
        // auto [barrier_f, barrier_g, barrier_H] = barrierEnergy_func.eval_with_derivatives(x); //
        
        // printf(" finding Barrier terms\n");
        // DenseMatrix<double> x_in_dense_format = vertex_data_to_matrix(tmp_geometry->inputVertexPositions);
        // auto [my_barrier_f, barrier_grad, barrier_hessian] = get_log_barrier_stuff(x_in_dense_format, Vector<double>::Ones(mesh->nVertices())); // thank u new c++
        // Eigen::SparseMatrix<double> my_barrier_H = tinyADify_barrier_hess(barrier_hessian);
        // Eigen::VectorX<double> my_barrier_g = tinyAD_flatten(barrier_grad);
        // TODO: smarter scheduling?
          
        // Eigen::SparseMatrix<double> total_H = bending_lambda * bending_H_proj + membrane_lambda * membrane_H_proj + CP_lambda * CP_H + barrier_lambda * barrier_H;
        // Eigen::VectorX<double> total_g = bending_lambda * bending_g + membrane_lambda * membrane_g + CP_lambda * CP_g + barrier_lambda * barrier_g;
        // double total_energy = bending_lambda * bending_f + membrane_lambda * membrane_f + CP_lambda * CP_energy + barrier_lambda * barrier_f;

        // TINYAD_INFO(" my barr vs TinyAD" << 
        //             "\n\t\t\tenergies: " << my_barrier_f << ", " << barrier_f <<
        //             "\n\t\t\tgrads: " << my_barrier_g.isApprox(barrier_g) << " "<< (barrier_g.transpose()-my_barrier_g.transpose()).norm()<<
        //             "\n\t\t\thessians:" << my_barrier_H.isApprox(barrier_H) << " "<< (barrier_H - my_barrier_H).norm());

        TINYAD_DEBUG_OUT("\t- Energy in iter " << i << ": bending= " << bending_lambda << "*" << bending_f << 
                         "\n\t\t\t\t membrane  = " << membrane_lambda << " * " << membrane_f <<
                         "\n\t\t\t\t CP        = " << CP_lambda << " * "<< CP_energy ); //<< 
                        //  "\n\t\t\t\t barr      = " << barrier_lambda << " * "<< barrier_f<< 
                        //  "\n\t\t\t\t\t total: " << total_energy);
        Eigen::VectorXd old_x = x;
        x = solve_QP_with_ineq((bending_lambda  * bending_H_proj + 
                                membrane_lambda * membrane_H_proj + 
                                reg_lambda      * reg_H +
                                CP_lambda       * CP_H)/2., //bending_lambda * bending_H_proj + membrane_lambda * membrane_H_proj + 
                                    bending_lambda  * (bending_g - bending_H_proj*x) + 
                                    membrane_lambda * (membrane_g - membrane_H_proj*x) + 
                                    reg_lambda      * (-2.*x) + 
                                    CP_lambda       * (CP_g - CP_H*x), old_x, //
                                constraint_matrix, constraint_rhs);
        // // Eigen::VectorXd d = TinyAD::newton_direction(total_g, total_H, solver);
        // shiftDiagonal(barrier_H, 1e-7); // TODO
        // printf("testing barrier H\n");
        // // LinearSolver<double>
        // try {
        //     PositiveDefiniteSolver<double> test_solver(barrier_H);
        // }
        // catch(const std::exception& e) {
        //     std::cerr << "barrier H NOT ok\n"<< e.what() << '\n';
        // }
        
        // shiftDiagonal(total_H, 1e-7); // TODO

        // Solver<double> newton_solver1(total_H);
        // Eigen::VectorXd d1 = newton_solver1.solve(-total_g);
        
        // PositiveDefiniteSolver<double> newton_solver(total_H);
        // Eigen::VectorXd d = newton_solver.solve(-total_g);
        // TINYAD_INFO(" resid: "<< (total_H * d + total_g).norm() << 
        //             "\n resid1: "<< (total_H * d1 + total_g).norm() <<
        //             "\n diff: "<< (d1 - d).norm());

        // x = TinyAD::line_search(x, d, total_energy, total_g,
        //                         [&] (const Eigen::VectorXd curr_x, bool print = false) {
        //                             if (!check_feasibility(unflat_tinyAD(curr_x)))
        //                                 return std::numeric_limits<double>::infinity();
        //                             double bending_e = bendingEnergy_func.eval(curr_x),
        //                                    membrane_e = membraneEnergy_func.eval(curr_x),
        //                                    CP_e = closest_point_energy(curr_x),
        //                                    log_barr_e = get_log_barrier_energy(unflat_tinyAD(curr_x));
        //                             double f_new = bending_lambda * bending_e + membrane_lambda * membrane_e + CP_lambda * CP_e + barrier_lambda * log_barr_e;  
        //                             if (print)
        //                                 TINYAD_WARNING("$$$$ Energies:"<<
        //                                                 "\n\t\t\t\t bending= " << bending_lambda << "*" << bending_e << 
        //                                                 "\n\t\t\t\tmembrane  = " << membrane_lambda << membrane_e <<
        //                                                 "\n\t\t\t\tCP  = " << CP_lambda << " * "<< CP_e<< 
        //                                                 "\n\t\t\t\tbarr= " << barrier_lambda << " * "<< log_barr_e<< 
        //                                                 "\n\t\t\t\t\t total:" << f_new);
        //                             return f_new;
        //                             },
        //                         0.5, 0.9, 200, 0.); // no clue whats good here for armijo constant; proly nothing since non-linear stuff happening
        //                       //smax, decay, max_iter, armijo constant
        
        // update tmp geometry
        bendingEnergy_func.x_to_data(x, [&] (Vertex v, const Eigen::Vector3d& p) {
            tmp_geometry->inputVertexPositions[v] = to_geometrycentral(p);
        });
        
        double step_norm = (x - old_x).norm();
        printf(" step norm is %9f\n", step_norm);

        // // dynamic remesh 
        if (dynamic_remesh){
            // joint_remesh(mesh, old_geometry, tmp_geometry, initial_mean_edge_len);
            split_only_remesh(mesh, old_geometry, tmp_geometry, initial_mean_edge_len);

            barrierEnergy_func = get_tinyAD_barrier_function();
            update_bending_rest_constants();
            bendingEnergy_func = get_tinyAD_bending_function();
            update_membrane_rest_constants();
            membraneEnergy_func = get_tinyAD_membrane_function();
            x = bendingEnergy_func.x_from_data([&] (Vertex v) {
                return to_eigen(tmp_geometry->inputVertexPositions[v.getIndex()]);
            });
        }
        
        // scheduled weights
        // if (step_norm < 0.05 && CP_energy > 0.05){
        //     printf(" Cp increment\n");    
        //     CP_lambda *= CP_mu;
        // }
        // else if (CP_energy < 0.05){
        //     printf(" barrier decrement\n");
        //     membrane_lambda *= CP_mu;
        //     bending_lambda *= CP_mu;
        //     barrier_lambda *= barrier_decay;
        // }
        barrier_lambda = get_scheduled_weight(init_barrier_lambda, final_barrier_lambda);
        bending_lambda = get_scheduled_weight(init_bending_lambda, final_bending_lambda);
        membrane_lambda = get_scheduled_weight(init_membrane_lambda, final_membrane_lambda);
        CP_lambda = get_scheduled_weight(init_CP_lambda, final_CP_lambda);
        // if (step_norm < convergence_eps)
        //     break;
    }
    DenseMatrix<double> new_points_mat = unflat_tinyAD(x); 
    return new_points_mat;
}



// // 
//  // joint remesh attempt
// joint_remesh(mesh, old_geometry, tmp_geometry, initial_mean_edge_len);
// // // update stuff after remesh; 
// // // ****** TODO how do I move this tinyAD stuff into a function??? ******
// printf(" re-compute old geometry constants\n");
// old_geometry->refreshQuantities();
// old_geometry->requireEdgeDihedralAngles();
// rest_dihedral_angles = old_geometry->edgeDihedralAngles;
// old_geometry->unrequireEdgeDihedralAngles();
// rest_bending_constant = get_bending_rest_constants();

// // re-compute bending function
// printf(" re-compute bending function\n");
// bendingEnergy_func = TinyAD::scalar_function<3>(mesh->vertices()); // 
// bendingEnergy_func.add_elements<4>(mesh->edges(), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
// {
//     using T = TINYAD_SCALAR_TYPE(element);

//     Edge e = element.handle;
    
//     if (e.isBoundary()) return (T)0.0;

//     Eigen::Vector3<T> p1 = element.variables(e.firstVertex());
//     Eigen::Vector3<T> p2 = element.variables(e.secondVertex());
    
//     Eigen::Vector3<T> p3 = element.variables(e.halfedge().next().tipVertex());
//     Eigen::Vector3<T> p4 = element.variables(e.halfedge().twin().next().tipVertex());

//     Eigen::Vector3<T> N1 = (p2 - p1).cross(p3 - p1); //faceNormals[e.halfedge().face()];
//     Eigen::Vector3<T> N2 = (p1 - p2).cross(p4 - p2);
//     Eigen::Vector3<T> edgeDir = (p2 - p1).normalized();
//     T dihedral_angle = atan2(edgeDir.dot(N1.cross(N2)), N1.dot(N2));
//     return (dihedral_angle - rest_dihedral_angles[e]) * (dihedral_angle - rest_dihedral_angles[e]) * rest_bending_constant[e];
// });
// // re-build x; since the size has changed 
// // tmp geo is updated in remeshing automatically
// x = bendingEnergy_func.x_from_data([&] (Vertex v) {
//     return to_eigen(tmp_geometry->inputVertexPositions[v.getIndex()]);
// });


///////// my old line search

// line search
// // printf(" line search\n");
// double _s_max = 0.5; // Initial step size
// double _shrink = 0.8;
// int _max_iters = 200; // 64
// double _armijo_const = 0.; //1e-4;
// bool try_one = _s_max > 1.0;

// Eigen::VectorX<double> x_new = x;
// double s = _s_max;
// double f_new = total_energy;
// double d_dot_g = d.dot(total_g);
// int j = 0;
// for (j = 0; j < _max_iters; ++j)
// {
//     x_new = x + s * d;
//     f_new = bendingEnergy_func.eval(x_new) + 
//                             CP_lambda * closest_point_energy(x_new);
//     if (f_new <= total_energy + _armijo_const * s * d_dot_g)
//         break; //  x new is good
//     else
//         s *= _shrink;
// }
// TINYAD_DEBUG_OUT("line ended at iter _" << j << "_ s: " << s << " \n                  fnew: " << f_new);
// x = x_new;