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
 #include "deformation.h"


// vector stuff
Eigen::Vector3d to_eigen(const geometrycentral::Vector3& _v){
    return Eigen::Vector3d(_v.x, _v.y, _v.z);
}

geometrycentral::Vector3 to_geometrycentral(
        const Eigen::Vector3d& _v) {
    return geometrycentral::Vector3 { _v.x(), _v.y() };
}

Vector<double> tinyAD_flatten(DenseMatrix<double> mat){
    size_t n = mat.rows();
    Vector<double> ans(n*3);
    for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j < n; j++)
            ans(3*i + j) = mat(i,j);
    }
    return ans;
}

Vector<double> vec32vec(Vector3 v){
    Vector<double> ans(3);
    ans[0] = v.x;
    ans[1] = v.y;
    ans[2] = v.z;
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
    closest_point_flat_operator = Eigen::SparseMatrix<double>(convex_mesh->nVertices(), 3*mesh->nVertices());
    std::vector<Eigen::Triplet<double>> tripletList, flat_tripletList;
    tripletList.reserve(3*mesh->nVertices());
    flat_tripletList.reserve(9*mesh->nVertices());

    for (Vertex c_v: convex_mesh->vertices()){
        Face closest_face;
        Edge closest_edge;
        Vector3 on_edge_projection;
        Vertex closest_vertex;
        double min_face_dist  = std::numeric_limits<double>::infinity(),
               min_edge_dist   = std::numeric_limits<double>::infinity(),
               min_vertex_dist  = std::numeric_limits<double>::infinity();

        Vector3 c_p = convex_geometry->inputVertexPositions[c_v];

        // there should be a smarter way of checking all elements in one loop; maybe even without a flagging them??
        for (Face f: mesh->faces()){ // check if face-projectable; while keeping track of closest vertex 
            Vector3 f_normal = new_geometry->faceNormal(f);
            double face_dist = abs(dot(f_normal, c_p - new_geometry->inputVertexPositions[f.halfedge().vertex()])); // using some point of f
            if (face_dist < min_face_dist){
                min_face_dist = face_dist;
                closest_face = f;
            }
        }
        for (Edge e: mesh->edges()){
            Vector3 A = convex_geometry->inputVertexPositions[e.firstVertex()], 
                    B = convex_geometry->inputVertexPositions[e.secondVertex()];
            Vector3 PB = B - c_p,
                    PA = A - c_p,
                    AB = B - A;
            Vector3 ortho_p = PB - AB*dot(AB, PB)/dot(AB,AB);
            double edge_dist =  ortho_p.norm();
            if (edge_dist < min_edge_dist){
                closest_edge = e;
                min_edge_dist = edge_dist;
                on_edge_projection = c_p + ortho_p;
            }
        }
        for (Vertex v: mesh->vertices()){
            Vector3 A = convex_geometry->inputVertexPositions[v];
            double vertex_dist = (A - c_p).norm();
            if (vertex_dist < min_vertex_dist){
                closest_vertex = v;
                min_vertex_dist = vertex_dist;
            }
        }

        // assign SurfacePoint and assign barycentric coordinates
        double min_dist = std::min(min_face_dist, std::min(min_edge_dist, min_vertex_dist));
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
            closest_point_assignment[c_v] = SurfacePoint(closest_face, bary_coor);

            // operator entries
            tripletList.emplace_back(c_v.getIndex(), v1.getIndex(), bary_coor.x);
            tripletList.emplace_back(c_v.getIndex(), v2.getIndex(), bary_coor.y);
            tripletList.emplace_back(c_v.getIndex(), v3.getIndex(), bary_coor.z);
            for (size_t d = 0; d < 3; d++){
                flat_tripletList.emplace_back(c_v.getIndex(), 3*v1.getIndex() + d, bary_coor.x);
                flat_tripletList.emplace_back(c_v.getIndex(), 3*v1.getIndex() + d, bary_coor.y);
                flat_tripletList.emplace_back(c_v.getIndex(), 3*v1.getIndex() + d, bary_coor.z);
            }
        }
        else if (min_dist == min_edge_dist){
            Vertex v1 = closest_edge.firstVertex(),
                   v2 = closest_edge.secondVertex();
            Vector3 A = convex_geometry->inputVertexPositions[v1], 
                    B = convex_geometry->inputVertexPositions[v2];
            double tVal = (on_edge_projection - A).norm()/(A - B).norm();
            closest_point_assignment[c_v] = SurfacePoint(closest_edge, tVal);

            // operator entries
            tripletList.emplace_back(c_v.getIndex(), v1.getIndex(), 1. - tVal);
            tripletList.emplace_back(c_v.getIndex(), v2.getIndex(), tVal);
            for (size_t d = 0; d < 3; d++){
                flat_tripletList.emplace_back(c_v.getIndex(), 3*v1.getIndex() + d, 1. - tVal);
                flat_tripletList.emplace_back(c_v.getIndex(), 3*v2.getIndex() + d, tVal);
            }
        }
        else {// min_dist == min_vertex_dist
            closest_point_assignment[c_v] = SurfacePoint(closest_vertex);   
            tripletList.emplace_back(c_v.getIndex(), closest_vertex.getIndex(), 1.);
            for (size_t d = 0; d < 3; d++)
                flat_tripletList.emplace_back(c_v.getIndex(), 3*closest_vertex.getIndex() + d, 1.);
        }
    }
    closest_point_operator.setFromTriplets(tripletList.begin(), tripletList.end());
    closest_point_flat_operator.setFromTriplets(flat_tripletList.begin(), flat_tripletList.end());
}


void DeformationSolver::build_constraint_matrix_and_rhs(){
    size_t nf = convex_mesh->nFaces();
    DenseMatrix<double> convex_face_normals(nf, 3);
    convex_geometry->requireFaceNormals(); // assuming outward normals
    constraint_matrix = face_data_to_matrix(convex_geometry->faceNormals);
    constraint_rhs = Vector<double>::Zero(nf);
}


VertexData<Vector3> DeformationSolver::solve_for_bending(VertexPositionGeometry* new_geometry){
    size_t num_var = 3 * mesh->nVertices();
    
    // Pre-compute constants from the old geometry; old \theta, e, h_e
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

    // Objective function
    auto bendingEnergy_func = TinyAD::scalar_function<3>(mesh->vertices());
    // Add objective term per edge
    bendingEnergy_func.add_elements<3>(mesh->edges(), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {
        // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);

        // Get variable 3D vertex positions
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

    // Assemble inital x vector from parametrization property.
    // x_from_data(...) takes a lambda function that maps
    // each variable handle to its initial 2D value (Eigen::Vector2d).
    int n = mesh->nVertices();
    Eigen::VectorXd x = bendingEnergy_func.x_from_data([&] (Vertex v) {
        return to_eigen(old_geometry->inputVertexPositions[v]);
    });

    // Projected Newton
    VertexPositionGeometry *tmp_geometry = new VertexPositionGeometry(*mesh);
    TinyAD::LinearSolver solver;
    int max_iters = 10;
    double convergence_eps = 1e-2,
           CP_lambda = 0.5;
    for (int i = 0; i < max_iters; ++i)
    {
        auto [f, g, H_proj] = bendingEnergy_func.eval_with_hessian_proj(x);

        bendingEnergy_func.x_to_data(x, [&] (Vertex v, const Eigen::Vector3d& p) {
            tmp_geometry->inputVertexPositions[v] = to_geometrycentral(p);
        });
        // directly get Closest Point energy stuff
        assign_closest_points(tmp_geometry);

        Eigen::SparseMatrix<double> A_CP = CP_lambda * closest_point_flat_operator;
        Eigen::SparseMatrix<double> new_H = H_proj + A_CP.transpose() * A_CP;
        Eigen::VectorX<double> CP_flat_grad = 2.*CP_lambda*tinyAD_flatten(closest_point_energy_gradient(tmp_geometry));
        Eigen::VectorX<double> new_g = g + CP_flat_grad;

        TINYAD_DEBUG_OUT("Energy in iteration " << i << ": " << f);
        Eigen::VectorXd d = TinyAD::newton_direction(g, new_H, solver);
        if (TinyAD::newton_decrement(d, g) < convergence_eps)
            break;

        // x = TinyAD::line_search(x, d, f, g, bendingEnergy_func);
        // manual line search for now
        double _s_max = 1.0; // Initial step size
        double _shrink = 0.8;
        int _max_iters = 64;
        double _armijo_const = 1e-4;
        bool try_one = _s_max > 1.0;

        Eigen::VectorX x_new;// = x;
        double s = _s_max;
        for (int i = 0; i < _max_iters; ++i)
        {
            x_new = _x0 + s * _d;
            const double f_new = _eval(x_new);
            TINYAD_ASSERT_EQ(f_new, f_new);
            if (armijo_condition(_f, f_new, s, _d, _g, _armijo_const))
                return x_new;

            if (try_one && s > 1.0 && s * _shrink < 1.0)
                s = 1.0;
            else
                s *= _shrink;
        }
    }
    // TINYAD_DEBUG_OUT("Final energy: " << func.eval(x));

    // // Write final x vector to parametrization property.
    // // x_to_data(...) takes a lambda function that writes the final value
    // // of each variable (Eigen::Vector2d) back our parametrization property.
    // func.x_to_data(x, [&] (Vertex v, const Eigen::Vector2d& p) {
    //     param[v] = to_geometrycentral(p);
    // });
    VertexData<Vector3> new_points(*mesh, Vector3::zero());
    // for (Vertex v : mesh->vertices()){
    //     new_points[v].x = X[v.getIndex()].get(GRB_DoubleAttr_X);
    //     new_points[v].y = X[v.getIndex() + nv].get(GRB_DoubleAttr_X);
    //     new_points[v].z = X[v.getIndex() + 2 * nv].get(GRB_DoubleAttr_X);
    // } 
    return new_points;
}