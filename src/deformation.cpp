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
    assign_closest_points(new_geometry); // make more efficient by not calling this all the time
    double energy = 0.;
    for (Vertex c_v: convex_mesh->vertices()){
        Vector3 c_p = convex_geometry->inputVertexPositions[c_v],
                closest_point = closest_point_assignment[c_v].interpolate(new_geometry->inputVertexPositions);
        energy += (c_p - closest_point).norm();
    }
    return energy;
}


VertexData<Vector3> DeformationSolver::closest_point_energy_gradient(VertexPositionGeometry *new_geometry){
    VertexData<Vector3> CP_energy_gradients(*mesh, Vector3::zero());

}


void DeformationSolver::assign_closest_points(VertexPositionGeometry *new_geometry){
    closest_point_assignment = VertexData<SurfacePoint>(*convex_mesh, SurfacePoint());
    
    // linear operator corresponding to assignments
    closest_point_operator = Eigen::SparseMatrix<double>(convex_mesh->nVertices(), mesh->nVertices());
    std::vector<Eigen::Triplet<double>> tripletList;

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
        }
        else {// min_dist == min_vertex_dist
            closest_point_assignment[c_v] = SurfacePoint(closest_vertex);   
            tripletList.emplace_back(c_v.getIndex(), closest_vertex.getIndex(), 1.);
        }
    }

    closest_point_operator.setFromTriplets(tripletList.begin(), tripletList.end());
}

