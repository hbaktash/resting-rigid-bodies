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
#include "geometrycentral/surface/meshio.h"


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
    // polyscope::show();
}


DeformationSolver::DeformationSolver(ManifoldSurfaceMesh *_mesh, VertexPositionGeometry *_old_geometry,
                                     ManifoldSurfaceMesh *_convex_mesh, VertexPositionGeometry *_convex_geometry,
                                     Vector3 _goal_G){
    mesh = _mesh;
    old_geometry = _old_geometry;
    convex_mesh = _convex_mesh;
    convex_geometry = _convex_geometry;
    goal_G = _goal_G;
}


DeformationSolver::DeformationSolver(ManifoldSurfaceMesh *_mesh, 
                                     VertexPositionGeometry *_old_geometry, VertexPositionGeometry *_deformed_geometry, 
                                     ManifoldSurfaceMesh *_convex_mesh, VertexPositionGeometry *_convex_geometry,
                                     Vector3 _goal_G){
    mesh = _mesh;
    old_geometry = _old_geometry;
    deformed_geometry = _deformed_geometry;
    convex_mesh = _convex_mesh;
    convex_geometry = _convex_geometry;
    goal_G = _goal_G;
}


double DeformationSolver::closest_point_energy(VertexPositionGeometry *new_geometry){
    double energy = 0.;
    DenseMatrix<double> new_pos_mat = vertex_data_to_matrix(new_geometry->inputVertexPositions),
                        convex_points_mat = vertex_data_to_matrix(convex_geometry->inputVertexPositions);
    DenseMatrix<double> closest_assignments = closest_point_operator*new_pos_mat;
    energy = (closest_assignments - convex_points_mat).squaredNorm();
    return energy;
}


double DeformationSolver::closest_point_energy(Vector<double> flat_new_pos_mat){
    DenseMatrix<double> convex_points_mat = vertex_data_to_matrix(convex_geometry->inputVertexPositions);
    return (closest_point_flat_operator*flat_new_pos_mat - tinyAD_flatten(convex_points_mat)).squaredNorm();
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
            // double call to get face->edge->vertex mappings
            closest_point_assignment[c_v] = get_robust_barycentric_point(closest_point_assignment[c_v], split_robustness_threshold); // 
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
            if (std::max({d1,d2,d3}) <= refinement_CP_threshold) {// or max?
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
            if (std::max({d1,d2}) <= refinement_CP_threshold){ // close and not robust
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
    Eigen::VectorX<bool> frozen_flags(3*mesh->nVertices());
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
        // vertex_only_assignment[cv] = true; // TODO make this an option
        split_occured = true;
        SurfacePoint assigned_cp = closest_point_assignment[cv];
        SurfacePoint robust_sp = get_robust_barycentric_point(assigned_cp, split_robustness_threshold); // hitting twice with robustness; could go face->vertex
        if (robust_sp.type != assigned_cp.type) // this doesn't happen anymore
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
    old_geometry->refreshQuantities();
    rest_dihedral_angles = old_geometry->edgeDihedralAngles;
    old_geometry->unrequireEdgeDihedralAngles();

    // rest e/h_e
    rest_bending_constant = EdgeData<double>(*mesh, 0.); // e/h_e

    for (Edge e : mesh->edges()) {
        if (e.isBoundary()) continue; // set to zero 
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
    old_geometry->refreshQuantities();
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

        // trying deformed constant rather than rest constant
        T current_edge_len_sqrd = (p2 - p1).squaredNorm();
        T current_area_sum = 0.5 * (N1.norm() + N2.norm());
        T current_bending_constant = current_edge_len_sqrd / current_area_sum;
        
        // rest area measure
        return (dihedral_angle - rest_dihedral_angles[e]) * (dihedral_angle - rest_dihedral_angles[e]) * rest_bending_constant[e];
        // deformed area measure
        // return (dihedral_angle - rest_dihedral_angles[e]) * (dihedral_angle - rest_dihedral_angles[e]) * current_bending_constant;
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

        // trying different W(I_rest_inv * I_new) 
        // -- conformal from MIPS: An Efficient Global Parametrization Method - Kai Hormann and G¨unther Greiner
        return rest_face_areas[f] * ((simil_matrix.transpose() * simil_matrix).trace() * 0.5/simil_matrix.determinant() -1.);
        // -- from "Discrete Shells" - Grinspun et al.
        // return rest_face_areas[f] * 
        //         (membrane_mu * simil_matrix.trace()/2. + membrane_lambda * simil_matrix.determinant()/4.
        //          -(membrane_mu/4. + membrane_lambda/2.) * log(simil_matrix.determinant()));
        // return rest_face_areas[f] * log((simil_matrix.transpose() * simil_matrix).trace()*0.5 / simil_matrix.determinant()); // 
        // return rest_face_areas[f] * current_triangle_area  * log((simil_matrix.transpose() * simil_matrix).trace()*0.5 / simil_matrix.determinant()); // 
        // return rest_face_areas[f] * current_triangle_area * (simil_matrix.transpose() * simil_matrix).trace()*0.5 / simil_matrix.determinant();
    });
    return membraneEnergy_func;
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


DenseMatrix<double> DeformationSolver::solve_for_bending(int visual_per_step, bool save_steps){
    size_t num_var = 3 * mesh->nVertices();
    std::cout << "num_var: " << num_var << std::endl;
    old_geometry->refreshQuantities();
    old_geometry->requireEdgeLengths();
    double initial_mean_edge_len = old_geometry->edgeLengths.toVector().mean();
    old_geometry->unrequireEdgeLengths();
    std::cout << "initial mean edge length: " << initial_mean_edge_len << std::endl;

    Eigen::SparseMatrix<double> cv_gauss_curve_diag;

    // tinyAD stuff
    std::cout << " initializing variables1\n";
    build_constraint_matrix_and_rhs();
    std::cout << " initializing variables2\n";
    update_bending_rest_constants();
    std::cout << " initializing variables3\n";
    auto bendingEnergy_func = get_tinyAD_bending_function();
    std::cout << " initializing variables4\n";
    update_membrane_rest_constants();
    std::cout << " initializing variables5\n";
    auto membraneEnergy_func = get_tinyAD_membrane_function();    
    std::cout << " initializing variables6\n";
    
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
        x *= 0.99;
    }

    VertexPositionGeometry *tmp_geometry = new VertexPositionGeometry(*mesh, unflat_tinyAD(x));
    
    // some parameters
    double barrier_lambda = init_barrier_lambda,
           bending_lambda = init_bending_lambda,
           membrane_lambda = init_membrane_lambda,
           CP_lambda = init_CP_lambda,
           G_lambda = 0;
    internal_pt = 1.;

    vertex_only_assignment = VertexData<bool>(*convex_mesh, false);
    frozen_assignment = VertexData<Vertex>(*convex_mesh, Vertex());

    double initial_LS_step_size = 1.;

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
        
        // update G and volume
        auto GV_pair = find_center_of_mass(*mesh, *tmp_geometry);
        current_G = GV_pair.first;
        current_volume = GV_pair.second;

        if (visual_per_step != 0 && i % visual_per_step == 0){
            printf(" visualizing step %d,  --mesh size: %d\n", i, mesh->nVertices());
            auto tmp_PSmesh = polyscope::registerSurfaceMesh("temp sol", tmp_geometry->inputVertexPositions, mesh->getFaceVertexList());
            // printf(" tmp geo container size %d\n", tmp_geometry->inputVertexPositions.size());
            // printf(" ps mesh size %d\n", tmp_PSmesh->vertexDataSize);
            polyscope::registerPointCloud("G curr tmp", std::vector<Vector3>({current_G}))->setPointColor({1.,0.,0.})->setPointRadius(0.01)->setEnabled(true);
            polyscope::registerPointCloud("G goal", std::vector<Vector3>({goal_G}))->setPointColor({0.,0.,1.})->setPointRadius(0.01)->setEnabled(true);
                
            tmp_PSmesh->setSurfaceColor({136./255., 229./255., 107./255.});
            tmp_PSmesh->setEdgeWidth(1.);
            tmp_PSmesh->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
            tmp_PSmesh->setEnabled(true);
            polyscope::screenshot();
            polyscope::frameTick();
            // polyscope::show();
        }
        if (save_steps){
            std::string output_name = "iter_" + std::to_string(i);
            std::string parent_path = "../meshes/hulls/nonconvex_deformation/deformation_saves/opt_seq/";
            std::string output_path = parent_path + std::string(output_name) + ".obj";
            // open file if doesnt exist
            if (!std::filesystem::exists(parent_path)){
                std::filesystem::create_directories(parent_path);
            }
            // write mesh
            std::cout << ANSI_FG_GREEN << " saving deformed mesh to " << output_path << ANSI_RESET << std::endl;
            writeSurfaceMesh(*mesh, *tmp_geometry, output_path);
        }

        // dynamic remesh 
        bool topo_changed = false;
        if (dynamic_remesh){
            // joint_remesh(mesh, old_geometry, tmp_geometry, initial_mean_edge_len);
            topo_changed = topo_changed || split_only_remesh(mesh, old_geometry, tmp_geometry, initial_mean_edge_len);
        }
        // split for small CP
        if (threshold_exceeded > 0){
            std::cout << ANSI_FG_RED << " **** splitting for small CP ****" << ANSI_RESET << std::endl;
            // modifies mesh, old_geometry, tmp_geometry
            topo_changed = topo_changed || split_barycentric_closest_points(tmp_geometry); 
        }
        // split when close to convex hull
        if (topo_changed){
            // update tineAD stuff
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
        double CP_energy = closest_point_energy(x);
        Eigen::VectorXd x_cv_flat = tinyAD_flatten(vertex_data_to_matrix(convex_geometry->inputVertexPositions));
        Eigen::SparseMatrix<double> A_CP = closest_point_flat_operator;
        Eigen::VectorX<double> CP_g = 2. * A_CP.transpose() * (A_CP * x - x_cv_flat); 
        Eigen::SparseMatrix<double> CP_H = 2. * A_CP.transpose() * A_CP;
        
        // get G diff energy and derivative
        // double Gdiff_f = (current_G - goal_G).norm2();
        // Eigen::VectorXd Gdiff_g = tinyAD_flatten(per_vertex_G_derivative(tmp_geometry));
        // Eigen::SparseMatrix<double> Gdiff_H = identityMatrix<double>(x.size());

        // x - x_0 regularizer
        // Eigen::VectorXd reg_g = Eigen::VectorXd::Zero(x.size());
        // Eigen::SparseMatrix<double> reg_H = 2.* identityMatrix<double>(x.size());

        // elastic stuff
        auto [bending_f, bending_g, bending_H_proj] = bendingEnergy_func.eval_with_hessian_proj(x); //
        auto [membrane_f, membrane_g, membrane_H_proj] = membraneEnergy_func.eval_with_hessian_proj(x); //
          
        Eigen::SparseMatrix<double> total_H = bending_lambda * bending_H_proj + membrane_lambda * membrane_H_proj + CP_lambda * CP_H;// + reg_lambda * reg_H + G_lambda * Gdiff_H; // + barrier_lambda * barrier_H;
        Eigen::VectorXd total_g = bending_lambda * bending_g + membrane_lambda * membrane_g + CP_lambda * CP_g ; // + G_lambda * Gdiff_g; // reg g is zero at x0 // + barrier_lambda * barrier_g;
        double total_energy = bending_lambda * bending_f + membrane_lambda * membrane_f + CP_lambda * CP_energy; // + G_lambda * Gdiff_f; // reg e is zero at x0 // barrier_lambda * barrier_f;
        

        std::cout << ANSI_FG_MAGENTA << "\t- Energy in iter " << i << ": bending= " << bending_lambda << "*" << bending_f << 
                         "\n\t\t\t\t membrane  = " << membrane_lambda << " * " << membrane_f <<
                        //  "\n\t\t\t\t G-diff    = " << G_lambda << " * " << Gdiff_f <<
                         "\n\t\t\t\t CP        = " << CP_lambda << " * "<< CP_energy << ANSI_RESET << std::endl << 
                         "\n\t\t\t\t\t total: " << bending_lambda * bending_f + membrane_lambda * membrane_f + CP_lambda * CP_energy << ANSI_RESET << std::endl;
        
        // polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("bending g", unflat_tinyAD(-bending_g))->setVectorColor({0.1,0.1,1})->setEnabled(false);
        // polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("membrane g", unflat_tinyAD(-membrane_g))->setVectorColor({1,0.1,0.1})->setEnabled(false);
        // polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("CP g", unflat_tinyAD(-CP_g))->setVectorColor({0.1,1,0.1})->setEnabled(false);
        // polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("total g wieghted", unflat_tinyAD(-total_g))->setVectorColor({0.5,0.5,0.5})->setEnabled(false);
        // polyscope::show();
        
        Eigen::VectorXd old_x = x;
        Eigen::VectorXd d;
        
        // use QP solver
        Eigen::MatrixX<bool> active_set = get_active_set_matrix(unflat_tinyAD(x), active_set_threshold);
        std::cout << ANSI_FG_YELLOW << " solving QP... " << ANSI_RESET << std::endl;
        // Eigen::VectorXd new_x = solve_QP_with_ineq_OSQP(total_H/2., total_g - total_H * x, old_x, //
        //                                             constraint_matrix, constraint_rhs, 
        //                                             get_frozen_flags(), get_frozen_x(),
        //                                             active_set);
        Eigen::VectorXd new_x = solve_QP_with_ineq_GRB(total_H/2., total_g - total_H * x, old_x, //
                                                    constraint_matrix, constraint_rhs, 
                                                    get_frozen_flags(), get_frozen_x(),
                                                    active_set);
        
        std::cout << ANSI_FG_GREEN << "  ... done solving QP" << ANSI_RESET << std::endl;
        d = new_x - old_x;
        
        // line search
        double opt_step = line_search(old_x, d, total_energy, total_g,
                        [&] (const Eigen::VectorXd curr_x, bool print = false) {
                                if (!check_feasibility(unflat_tinyAD(curr_x)))
                                    return std::numeric_limits<double>::infinity();
                                double bending_e = bendingEnergy_func.eval(curr_x),
                                    membrane_e = membraneEnergy_func.eval(curr_x),
                                    CP_e = get_raw_CP_energy(*mesh, curr_x, *convex_mesh, *convex_geometry), // assigns CP at every step
                                    reg_e = 0;  //(curr_x - old_x).squaredNorm();
                                double f_new = bending_lambda * bending_e + membrane_lambda * membrane_e + CP_lambda * CP_e + reg_lambda * reg_e; // + barrier_lambda * log_barr_e;  
                                return f_new;
                            },
                        initial_LS_step_size, 0.9, 200);
                    // smax, decay, max_iter
        x += opt_step * d;
        // initial_LS_step_size = opt_step != 0 ? 10. * opt_step : 1.; // adds 22 steps to the current step size
        std::cout << ANSI_FG_YELLOW << "step size: " << opt_step << " next init LS step " << initial_LS_step_size << ANSI_RESET << std::endl;
        // update tmp geometry
        VertexData<Vector3> old_positions = tmp_geometry->inputVertexPositions;
        bendingEnergy_func.x_to_data(x, [&] (Vertex v, const Eigen::Vector3d& p) {
            tmp_geometry->inputVertexPositions[v] = to_geometrycentral(p);
        });

        double smallest_dist = 1e9, largest_dist = 0.;
        for (Vertex v: tmp_geometry->mesh.vertices()){
            double dist = (tmp_geometry->inputVertexPositions[v] - old_positions[v]).norm();
        }
        double step_norm = (x - old_x).norm();
        std::cout << ANSI_FG_YELLOW << "step norm: " << step_norm << ANSI_RESET << std::endl;
        // internal_pt *= internal_growth_p;
        
        if (step_norm < 0.06){ // 0.06 for conformal
            if (CP_lambda == final_CP_lambda){ // && G_lambda == final_G_lambda
                break;
            }
            if (get_raw_CP_energy(*mesh, x, *convex_mesh, *convex_geometry) < 1e-5){
                break;
            }
            CP_lambda = CP_lambda > final_CP_lambda ? final_CP_lambda : CP_lambda/internal_growth_p;
            bending_lambda = bending_lambda > final_bending_lambda ? final_bending_lambda : bending_lambda/(std::min({2. * internal_growth_p, (1.+internal_growth_p)/2.}));
            // G_lambda = G_lambda > final_G_lambda ? final_G_lambda : G_lambda/internal_growth_p;
            std::cout << ANSI_FG_RED << " step norm is small, increasing CP lambda to " << CP_lambda << ANSI_RESET << std::endl;
        }
        // bending_lambda = get_scheduled_weight(init_bending_lambda, final_bending_lambda);
        // membrane_lambda = get_scheduled_weight(init_membrane_lambda, final_membrane_lambda);
        // barrier_lambda = get_scheduled_weight(init_barrier_lambda, final_barrier_lambda);
        // if (step_norm < convergence_eps)
        //     break;
        // if (step_norm == 0 && bending_lambda == final_bending_lambda && membrane_lambda == final_membrane_lambda && CP_lambda == final_CP_lambda && barrier_lambda == final_barrier_lambda)
        //     break;
    }
    // last tick
    auto tmp_PSmesh = polyscope::registerSurfaceMesh("temp sol", tmp_geometry->inputVertexPositions, mesh->getFaceVertexList());
    tmp_PSmesh->setSurfaceColor({136./255., 229./255., 107./255.});
    tmp_PSmesh->setEdgeWidth(1.);
    tmp_PSmesh->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
    tmp_PSmesh->setEnabled(true);
    // polyscope::screenshot();
    polyscope::frameTick();
    // save deformed shape
    

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
                                                   bool energy_plot, int* current_iter, float** energies_log,
                                                   bool save_steps){
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
    Eigen::VectorXd x0 = bendingEnergy_func.x_from_data([&] (Vertex v) {
        return to_eigen(old_geometry->inputVertexPositions[v.getIndex()]);
    });
    VertexPositionGeometry *tmp_geometry = new VertexPositionGeometry(*mesh, unflat_tinyAD(x));
    assert(check_feasibility(unflat_tinyAD(x))); // the deformed shape should be feasible
    std::cout << "temp geo vs def geo: \n -- "<< tmp_geometry->mesh.nVertices() << "\n -- " << deformed_geometry->mesh.nVertices() << "\n";

    // some parameters
    double bending_lambda = init_bending_lambda,
           membrane_lambda = init_membrane_lambda,
           G_lambda = init_G_lambda;
    internal_pt = 1.;

    // trying reg only

    // 
    Eigen::VectorXd flat_dist_mult = flat_distance_multiplier(tmp_geometry, true);
    flat_dist_mult /= flat_dist_mult.maxCoeff();
    for(int i = 0; i < flat_dist_mult.size(); ++i)
        flat_dist_mult[i] = sqrt(flat_dist_mult[i]);

    for (int i = 0; i < filling_max_iter; ++i) {    
        // update G
        auto GV_pair = find_center_of_mass(*mesh, *tmp_geometry);
        current_G = GV_pair.first;
        current_volume = GV_pair.second;
        // visuals
        if (visual_per_step != 0){
            if ((i+1) % visual_per_step == 0){
                printf(" visualizing step %d\n", i);
                polyscope::registerPointCloud("G current iter", std::vector<Vector3>({current_G}))->setPointColor({1.,0.,0.})->setPointRadius(0.01)->setEnabled(true);
                polyscope::registerPointCloud("G goal", std::vector<Vector3>({goal_G}))->setPointColor({0.,0.,1.})->setPointRadius(0.01)->setEnabled(true);
                
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
        if (save_steps){
            std::string output_name = "iter_" + std::to_string(i);
            std::string parent_path = "../meshes/hulls/nonconvex_deformation/deformation_saves/G_opt_seq/";
            std::string output_path = parent_path + std::string(output_name) + ".obj";
            // open file if doesnt exist
            if (!std::filesystem::exists(parent_path)){
                std::filesystem::create_directories(parent_path);
            }
            // write mesh
            std::cout << ANSI_FG_GREEN << " saving deformed mesh to " << output_path << ANSI_RESET << std::endl;
            writeSurfaceMesh(*mesh, *tmp_geometry, output_path);
        }

        // get G diff energy and derivative
        double Gdiff_f = (current_G - goal_G).norm2();
        Eigen::VectorXd Gdiff_g = tinyAD_flatten(per_vertex_G_derivative(tmp_geometry));

        // ||x - x_0|| regularizer
        double reg_f = (x - x0).squaredNorm();
        Eigen::VectorXd reg_g = 2.*(x - x0);
        reg_lambda = 0.;
        // Eigen::SparseMatrix<double> reg_H = 2.* identityMatrix<double>(x.size());
        
        // elastic stuff; comment out for peformance when their lambda is 0
        // auto [bending_f, bending_g] = bendingEnergy_func.eval_with_gradient(x); //
        // auto [membrane_f, membrane_g] = membraneEnergy_func.eval_with_gradient(x); //

        Eigen::VectorXd total_g = G_lambda * Gdiff_g + reg_lambda * reg_g; // + bending_lambda * bending_g + membrane_lambda * membrane_g; // reg g is zero at x0 // + barrier_lambda * barrier_g;
        double total_energy = G_lambda * Gdiff_f + reg_lambda * reg_f;     // + bending_lambda * bending_f + membrane_lambda * membrane_f; // reg e is zero at x0 // barrier_lambda * barrier_f;


        if (!use_static_dists){
            flat_dist_mult = flat_distance_multiplier(tmp_geometry, true);
            flat_dist_mult /= flat_dist_mult.maxCoeff();
            for(int i = 0; i < flat_dist_mult.size(); ++i)
                flat_dist_mult[i] = sqrt(flat_dist_mult[i]);
        }
        
        // DEBUG
        // printf("temp sol size %d\n other shit size", polyscope::getSurfaceMesh("temp sol")->vertexDataSize, unflat_tinyAD(flat_dist_mult).col(0).size());
        // polyscope::getSurfaceMesh("temp sol")->addVertexScalarQuantity("normalized distance to CV", unflat_tinyAD(flat_dist_mult).col(0))->setEnabled(false);
        // polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("total g un-wieghted", -1*unflat_tinyAD(total_g))->setEnabled(false);
        
        // // sobolev diffuse grads
        Eigen::VectorXd diffused_total_g = tinyAD_flatten(sobolev_diffuse_gradients(unflat_tinyAD(total_g), *mesh, *tmp_geometry, G_deform_sobolev_lambda, sobolev_p));
        
        total_g = diffused_total_g.cwiseProduct(flat_dist_mult.cwiseAbs2());
        // total_g = total_g.cwiseProduct(flat_dist_mult.cwiseAbs2());

        // std::cout << "temp sol size:" << mesh->nVertices() << " gghat size: " << G_Ghat_g.size() << "unflat gghat size" << unflat_tinyAD(G_Ghat_g).rows() << std::endl;
        // polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("total g diffused", -1*unflat_tinyAD(diffused_total_g))->setEnabled(false);
        // polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("bending g", -1*unflat_tinyAD(bending_g))->setEnabled(false);
        // polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("membrane g", -1*unflat_tinyAD(membrane_g))->setEnabled(false);        
        // polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("reg g unW", -1*unflat_tinyAD(reg_g))->setEnabled(false);
        // polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("Gdiff unW", -1*unflat_tinyAD(Gdiff_g))->setEnabled(false);
        // polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("reg g W", -1*unflat_tinyAD(reg_g.cwiseProduct(flat_dist_mult.cwiseAbs2())))->setEnabled(true);
        // polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("Gdiff w", -1*unflat_tinyAD(Gdiff_g.cwiseProduct(flat_dist_mult.cwiseAbs2())))->setEnabled(true);
        // polyscope::getSurfaceMesh("temp sol")->addVertexVectorQuantity("total g wieghted", -1*unflat_tinyAD(total_g))->setEnabled(true);

        // if (i % 100 == 1){
        //     polyscope::show();
        // }

        std::cout << ANSI_FG_MAGENTA << 
                         "\t- Energy in iter " << i << 
                        //  "\n\t\t\t\t bending= " << bending_lambda << "*" << bending_f << 
                        //  "\n\t\t\t\t membrane  = " << membrane_lambda << " * " << membrane_f <<
                         "\n\t\t\t\t G-diff    = " << G_lambda  << " * " << Gdiff_f <<
                         "\n\t\t\t\t reg       = " << reg_lambda << " * " << reg_f << ANSI_RESET << std::endl <<
                         "\n\t\t\t\t\t total:  " << // bending_lambda * bending_f + membrane_lambda * membrane_f + 
                                                   G_lambda * Gdiff_f + reg_lambda * reg_f << ANSI_RESET << std::endl;
        
        Eigen::VectorXd old_x = x;
        Eigen::VectorXd d = -total_g;
        if (use_QP_solver){
            Eigen::MatrixX<bool> active_set = get_active_set_matrix(unflat_tinyAD(x), active_set_threshold);
            Eigen::VectorX<bool> frozen_flags(x.size());
            frozen_flags.fill(false);

            Eigen::SparseMatrix<double> Q = identityMatrix<double>(num_var);
            Eigen::VectorXd c = -2. * (x + 1. * d); // -2x_0

            std::cout << ANSI_FG_YELLOW << " solving QP... " << ANSI_RESET << std::endl;
            Eigen::VectorXd new_x = solve_QP_with_ineq_GRB(Q, c, old_x, //
                                                        constraint_matrix, constraint_rhs, 
                                                        frozen_flags, old_x,
                                                        active_set);
            std::cout << ANSI_FG_GREEN << "  ... done solving QP" << ANSI_RESET << std::endl;
            d = new_x - old_x;
            // x = new_x;
        }
        std::cout << ANSI_FG_YELLOW << "d norm(): " << d.norm() << ANSI_RESET << std::endl;
        double initial_LS_step_size = 1;
        double opt_step = line_search(old_x, d, total_energy, total_g,
                        [&] (const Eigen::VectorXd curr_x, bool print = false) {
                                if (!check_feasibility(unflat_tinyAD(curr_x)))
                                    return std::numeric_limits<double>::infinity();
                                // find new G at every step
                                VertexPositionGeometry *tmp_tmp_geo = new VertexPositionGeometry(*mesh, unflat_tinyAD(curr_x));
                                Vector3 new_G = find_center_of_mass(*mesh, *tmp_tmp_geo).first;
                                
                                double energy = 0;
                                // energy += bending_lambda * bendingEnergy_func.eval(curr_x);
                                // energy += membrane_lambda * membraneEnergy_func.eval(curr_x);
                                energy += G_lambda * (new_G - goal_G).norm2();
                                energy += reg_lambda * (curr_x - x0).squaredNorm();
                                return energy;
                            },
                        initial_LS_step_size, 0.95, 200);
                    // smax, decay, max_iter, armijo constant
        x += opt_step * d;
        // update tmp geometry
        bendingEnergy_func.x_to_data(x, [&] (Vertex v, const Eigen::Vector3d& p) {
            tmp_geometry->inputVertexPositions[v] = to_geometrycentral(p);
        });

        std::cout << ANSI_FG_YELLOW << "step norm(): " << opt_step << ANSI_RESET << std::endl;
        if ((x - old_x).norm() <= 1e-6){
            G_lambda = G_lambda/internal_growth_p;
            if (G_lambda > final_G_lambda)
                G_lambda = final_G_lambda;
            std::cout << ANSI_FG_RED << " step norm is small, increasing G lambda to " << G_lambda << ANSI_RESET << std::endl;
        }
        if (Gdiff_f < 1e-4){ // converged
            break;
        }
        // internal_pt *= internal_growth_p;
        // G_lambda = get_scheduled_weight(init_G_lambda, final_G_lambda);
    }
    DenseMatrix<double> new_points_mat = unflat_tinyAD(x);
    deformed_geometry = new VertexPositionGeometry(*mesh, new_points_mat);
    return new_points_mat;
}
