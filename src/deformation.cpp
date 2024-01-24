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
#include "polyscope/surface_mesh.h"
#include "deformation.h"


// vector stuff
Eigen::Vector3d to_eigen(const geometrycentral::Vector3& _v){
    return Eigen::Vector3d(_v.x, _v.y, _v.z);
}

geometrycentral::Vector3 to_geometrycentral(
        const Eigen::Vector3d& _v) {
    return geometrycentral::Vector3 { _v.x(), _v.y(), _v.z() };
}

Vector<double> tinyAD_flatten(DenseMatrix<double> mat){
    size_t n = mat.rows();
    Vector<double> ans(n*3);
    for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j < 3; j++)
            ans(3*i + j) = mat(i,j);
    }
    return ans;
}

DenseMatrix<double> unflat_tinyAD(Vector<double> flat_mat){
    size_t n = flat_mat.size()/3;
    DenseMatrix<double> mat(n, 3);
    for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j < 3; j++)
            mat.coeffRef(i,j) = flat_mat(3*i + j);
    }
    return mat;
}

SparseMatrix<double> tinyADify_barrier_hess(std::vector<DenseMatrix<double>> hessians){
    size_t n = hessians.size();
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(9*n);

    for (size_t i = 0; i < n; i++){
        DenseMatrix<double> hess_i = hessians[i];
        for (size_t k = 0; k < 3; k++){
            for (size_t l = 0; l < 3; l++)
                tripletList.emplace_back(3 * i + k, 3 * i + l, hess_i(k, l));
        }
    }
    Eigen::SparseMatrix<double> tinyADfied_hess(3*n, 3*n);
    tinyADfied_hess.setFromTriplets(tripletList.begin(), tripletList.end());
    return tinyADfied_hess;
}


Vector<double> vec32vec(Vector3 v){
    Vector<double> ans(3);
    ans[0] = v.x;
    ans[1] = v.y;
    ans[2] = v.z;
    return ans;
}

Vector3 vec_to_GC_vec3(Vector<double> vec){
    Vector3 ans;
    ans.x = vec[0];
    ans.y = vec[1];
    ans.z = vec[2];
    return ans;
}

DenseMatrix<double> vertex_data_to_matrix(VertexData<Vector3> positions){
    size_t n = positions.getMesh()->nVertices();
    DenseMatrix<double> mat(n, 3);
    for (Vertex v: positions.getMesh()->vertices()){
        Vector3 p = positions[v];
        mat.row(v.getIndex()) = vec32vec(p);
    }
    return mat;
}


DenseMatrix<double> face_data_to_matrix(FaceData<Vector3> fdata){
    size_t n = fdata.getMesh()->nFaces();
    DenseMatrix<double> mat(n, 3);
    for (Face f: fdata.getMesh()->faces()){
        Vector3 p = fdata[f];
        mat.row(f.getIndex()) = vec32vec(p);
    }
    return mat;
}

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
    energy = (closest_point_operator*new_pos_mat - convex_points_mat).norm();
    // for (Vertex c_v: convex_mesh->vertices()){
    //     Vector3 c_p = convex_geometry->inputVertexPositions[c_v],
    //             closest_point = closest_point_assignment[c_v].interpolate(new_geometry->inputVertexPositions);
    //     energy += (c_p - closest_point).norm();
    // }
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


void DeformationSolver::assign_closest_points(VertexPositionGeometry *new_geometry){
    closest_point_assignment = VertexData<SurfacePoint>(*convex_mesh, SurfacePoint());
    
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
        double min_dist = std::min(min_face_dist, std::min(min_edge_dist, min_vertex_dist));
        if (min_dist == min_face_dist){ // assuming triangle mesh
            // printf("closest is a face\n");
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
            closest_point_assignment[c_v] = SurfacePoint(closest_face, bary_coor);

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
            Vertex v1 = closest_edge.firstVertex(),
                   v2 = closest_edge.secondVertex();
            Vector3 A = new_geometry->inputVertexPositions[v1], 
                    B = new_geometry->inputVertexPositions[v2];
            double tVal = (on_edge_projection - A).norm()/(A - B).norm();
            closest_point_assignment[c_v] = SurfacePoint(closest_edge, tVal);

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
            closest_point_assignment[c_v] = SurfacePoint(closest_vertex);   
            tripletList.emplace_back(c_v.getIndex(), closest_vertex.getIndex(), 1.);
            for (size_t d = 0; d < 3; d++)
                flat_tripletList.emplace_back(3 * c_v.getIndex() + d, 3*closest_vertex.getIndex() + d, 1.);
        }
    }
    closest_point_operator.setFromTriplets(tripletList.begin(), tripletList.end());
    closest_point_flat_operator.setFromTriplets(flat_tripletList.begin(), flat_tripletList.end());
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
DeformationSolver::get_log_barrier_stuff(DenseMatrix<double> new_pos_mat){
    double energy = 0.;
    size_t n = new_pos_mat.rows(),
           nf = convex_mesh->nFaces();
    assert(new_pos_mat.cols() == 3);
    DenseMatrix<double> NP = constraint_matrix * new_pos_mat.transpose(); // nf by n; col j is N*P_j
    DenseMatrix<double> repeated_rhs = constraint_rhs.replicate(1, n);
    assert(repeated_rhs.cols() == n); // same as the number of new pos vertices
    DenseMatrix<double> diff_ij = repeated_rhs - NP; // col j is N*P_j - rhs; i.e. -f_i(P_j) which should be positive
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
        hessians[j] = DenseMatrix<double>::Zero(3,3);
        for (size_t i = 0; i < nf; i++){
            DenseMatrix<double> g_gT = grads.row(j).transpose() * grads.row(j);
            if(!(g_gT.cols() == 3 && g_gT.rows() == 3))
                throw std::logic_error(" gt size "+ std::to_string(g_gT.rows()) + "," + std::to_string(g_gT.cols()) + "\n");
            double diff_ij_sqr = diff_ij.coeff(i,j) * diff_ij.coeff(i,j);
            hessians[j] += diff_ij_sqr * g_gT; // sum over N_i / (diff_ij)
        }
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


DenseMatrix<double> DeformationSolver::solve_for_bending(int visual_per_step){
    size_t num_var = 3 * mesh->nVertices();
    
    // Pre-compute constants from the old geometry; old \theta, e, h_e
    printf(" pre computing rest constants\n");
    old_geometry->requireEdgeDihedralAngles();
    EdgeData<double> rest_dihedral_angles = old_geometry->edgeDihedralAngles;
    old_geometry->unrequireEdgeDihedralAngles();
    EdgeData<double> rest_constant(*mesh); // e/h_e
    for (Edge e : mesh->edges()) {
        // Get 3D vertex positions
        Vector3 p1 = old_geometry->inputVertexPositions[e.firstVertex()];
        Vector3 p2 = old_geometry->inputVertexPositions[e.secondVertex()];
        // length and area
        double e_len_sqrd = (p1 - p2).norm2();
        double area1 = old_geometry->faceArea(e.halfedge().face()),
               area2 = old_geometry->faceArea(e.halfedge().twin().face());
        // Save 2-by-2 matrix with edge vectors as colums
        rest_constant[e] = e_len_sqrd/(area1+area2);
    };

    printf(" adding objective terms\n");
    // Objective function
    auto bendingEnergy_func = TinyAD::scalar_function<3>(mesh->vertices());
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

        return (dihedral_angle - rest_dihedral_angles[e]) * (dihedral_angle - rest_dihedral_angles[e]) * rest_constant[e];
    });


    printf(" initializing variables\n");
    int n = mesh->nVertices();
    Eigen::VectorXd x = bendingEnergy_func.x_from_data([&] (Vertex v) {
        return to_eigen(old_geometry->inputVertexPositions[v.getIndex()]);
    });

    build_constraint_matrix_and_rhs();
    //TODO: enforce feasibility
    while (!check_feasibility(unflat_tinyAD(x))){ // assuming centered; which is
        printf("  -- scaling for feasibility -- \n");
        x *= 0.9;
    }

    // Projected Newton
    VertexPositionGeometry *tmp_geometry = new VertexPositionGeometry(*mesh);
    TinyAD::LinearSolver solver;
    double convergence_eps = 1e-5;
    double barrier_lambda = barrier_init_lambda;
    printf(" starting opt\n");
    // polyscope::frameTick();
    auto tmp_PSmesh = polyscope::registerSurfaceMesh("temp sol", tmp_geometry->inputVertexPositions, mesh->getFaceVertexList());
    tmp_PSmesh->setSurfaceColor({136./255., 229./255., 107./255.});
    tmp_PSmesh->setEdgeWidth(1.);
    tmp_PSmesh->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
    tmp_PSmesh->setEnabled(true);
    // polyscope::show(); // close window to start animations
    for (int i = 0; i < filling_max_iter; ++i) {
        // visuals
        if (visual_per_step != 0){
            if (i % visual_per_step == 0){
                printf(" visualizing step %d\n", i);
                // auto tmp_PSmesh = polyscope::registerSurfaceMesh("temp sol "+std::to_string(i), tmp_geometry->inputVertexPositions, mesh->getFaceVertexList());
                // tmp_PSmesh->setSurfaceColor({136./255., 229./255., 107./255.});
                // tmp_PSmesh->setEdgeWidth(1.);
                // tmp_PSmesh->setBackFacePolicy(polyscope::BackFacePolicy::Custom);
                // tmp_PSmesh->setEnabled(false);
                tmp_PSmesh->updateVertexPositions(tmp_geometry->inputVertexPositions);
                polyscope::frameTick();
            }
        }

        printf(" Hessian eval\n ");
        auto [bending_f, bending_g, bending_H_proj] = bendingEnergy_func.eval_with_hessian_proj(x); //

        bendingEnergy_func.x_to_data(x, [&] (Vertex v, const Eigen::Vector3d& p) { // don't really need this
            tmp_geometry->inputVertexPositions[v] = to_geometrycentral(p);
        });
        // directly get Closest Point energy stuff
        if (!one_time_CP_assignment || i == 0){
            printf(" assigning  closest points\n");
            assign_closest_points(tmp_geometry);
        }
        printf(" finding CP terms\n");
        // CP stuff
        double CP_energy = closest_point_energy(tmp_geometry);
        Eigen::SparseMatrix<double> A_CP = closest_point_flat_operator;
        Eigen::VectorX<double> CP_flat_g = 2.*tinyAD_flatten(closest_point_energy_gradient(tmp_geometry));
        Eigen::SparseMatrix<double> H_CP = A_CP.transpose() * A_CP;
        
        // barrier stuff
        printf(" finding Barrier terms\n");
        DenseMatrix<double> x_in_dense_format = vertex_data_to_matrix(tmp_geometry->inputVertexPositions);
        auto [barrier_energy, barrier_grad, barrier_hessian] = get_log_barrier_stuff(x_in_dense_format); // thank u new c++
        Eigen::SparseMatrix<double> barrier_H = tinyADify_barrier_hess(barrier_hessian);
        Eigen::VectorX<double> barrier_g = tinyAD_flatten(barrier_grad);
        // TODO: what to do about scheduling
        // Given a frozen CP lambda at iter i; solve with constraints
            
        size_t inner_max_iters = 50;
        
        // TODO: make sure x is feasible
        // for (size_t inner_iter = 0; inner_iter < inner_max_iters; inner_iter++){
            
        Eigen::SparseMatrix<double> total_H = bending_H_proj + CP_lambda * H_CP + barrier_lambda * barrier_H;
        Eigen::VectorX<double> total_g = bending_g + CP_lambda * CP_flat_g + barrier_lambda * barrier_g;
        double total_energy = bending_f + CP_lambda * CP_energy + barrier_lambda * barrier_energy;

        TINYAD_DEBUG_OUT("\t- Energy in iter " << i << ": bending= " << bending_f << 
                            "\n\t\t\t\tCP  = " << CP_lambda << " * "<< CP_energy<< 
                            "\n\t\t\t\tbarr= " << barrier_lambda << " * "<< barrier_energy<< 
                            "\n\t\t\t\t\t total:" << total_energy);
        Eigen::VectorXd d = TinyAD::newton_direction(total_g, total_H, solver);
        // if (visual_per_step != 0){
        //     if (i % visual_per_step == 0){
        //         printf(" visualizing interior step %d\n", inner_iter);
        //           polyscope::registerSurfaceMesh("interior temp sol "+std::to_string(i), unflat_tinyAD(x), mesh->getFaceVertexList())->setEnabled(false);
        //     }
        // }
        Eigen::VectorXd old_x = x;
        x = TinyAD::line_search(x, d, total_energy, total_g,
                                [&] (const Eigen::VectorXd curr_x) {
                                    if (!check_feasibility(unflat_tinyAD(curr_x)))
                                        return std::numeric_limits<double>::infinity();
                                    return bendingEnergy_func.eval(curr_x) + 
                                            CP_lambda * closest_point_energy(curr_x) + 
                                            barrier_lambda * get_log_barrier_energy(unflat_tinyAD(curr_x));
                                    },
                                0.5, 0.8, 64, 0.); // no clue whats good here for armijo constant
                              //smax, decay, max_iter, armijo constant
        // }
        // scheduled weights
        barrier_lambda *= barrier_decay;
        CP_lambda *= CP_mu;
        double step_norm = (x - old_x).norm();
        printf(" step norm is %f\n", step_norm);
        if (step_norm < convergence_eps)
            break;
    }
    // polyscope::show();
    DenseMatrix<double> new_points_mat = unflat_tinyAD(x); 
    return new_points_mat;
}






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