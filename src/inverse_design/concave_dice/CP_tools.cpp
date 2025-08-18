#include "inverse_design/concave_dice/CP_tools.h"



double get_raw_CP_energy(ManifoldSurfaceMesh &concave_mesh, Eigen::VectorXd concave_poses_tiny_AD_flattened, 
                       ManifoldSurfaceMesh &convex_mesh, VertexPositionGeometry &convex_geometry){
    
    VertexPositionGeometry concave_geometry = VertexPositionGeometry(concave_mesh, unflat_tinyAD(concave_poses_tiny_AD_flattened));
    
    double cp_energy = 0.;

    for (Vertex c_v: convex_mesh.vertices()){
        Vector3 on_edge_projection;
        double min_face_dist  = std::numeric_limits<double>::infinity(),
               min_edge_dist   = std::numeric_limits<double>::infinity(),
               min_vertex_dist  = std::numeric_limits<double>::infinity();

        Vector3 c_p = convex_geometry.inputVertexPositions[c_v];

        // there should be a smarter way of checking all elements in one loop; maybe even without a flagging them??
        for (Face f: concave_mesh.faces()){ 
            Vector3 f_normal = concave_geometry.faceNormal(f);
            Halfedge curr_he  = f.halfedge(),
                    first_he = f.halfedge();
            // assume outward normals??
            double face_dist = dot(f_normal, c_p - concave_geometry.inputVertexPositions[f.halfedge().vertex()]); // using some point of f
            if (face_dist < 0) continue; // not on the right side of the face (outward normal)
            
            bool face_is_projectable = true;
            while (true){ // checking if face-projectable
                Vertex v1 = curr_he.tailVertex(), v2 = curr_he.tipVertex();
                Vector3 A = concave_geometry.inputVertexPositions[v1], 
                        B = concave_geometry.inputVertexPositions[v2];
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
            }
        }
        for (Edge e: concave_mesh.edges()){
            Vector3 A = concave_geometry.inputVertexPositions[e.firstVertex()], 
                    B = concave_geometry.inputVertexPositions[e.secondVertex()];
            Vector3 PB = B - c_p,
                    PA = A - c_p,
                    AB = B - A;
            bool edge_is_projectable = dot(PB, AB) >= 0 && dot(PA, -AB) >= 0;
            if (edge_is_projectable){
                Vector3 ortho_p = PB - AB*dot(AB, PB)/dot(AB,AB);
                double edge_dist =  ortho_p.norm();
                if (edge_dist < min_edge_dist){
                    min_edge_dist = edge_dist;
                    on_edge_projection = c_p + ortho_p;
                }
            }
        }
        for (Vertex v: concave_mesh.vertices()){
            Vector3 A = concave_geometry.inputVertexPositions[v];
            double vertex_dist = (A - c_p).norm();
            if (vertex_dist < min_vertex_dist){
                min_vertex_dist = vertex_dist;
            }
        }

        // assign SurfacePoint and assign barycentric coordinates
        // printf(" at cv %d fd %f, ed %f, vd %f\n", c_v.getIndex(), min_face_dist, min_edge_dist, min_vertex_dist);
        double min_dist = std::min({min_vertex_dist, min_edge_dist, min_face_dist}); // ordering is important in case of equality
        cp_energy += min_dist * min_dist;
    }

    return cp_energy;
}