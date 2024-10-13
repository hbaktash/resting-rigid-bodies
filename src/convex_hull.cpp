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
#include "convex_hull.h"


std::tuple<std::vector<std::vector<size_t>>, std::vector<size_t>, std::vector<Vector3>>
get_convex_hull(Eigen::MatrixX3d point_cloud){
    size_t n = point_cloud.rows();
    size_t dim = 3;
    std::vector<Vector3> point_set;
    for (size_t i = 0; i < n; i++){
        Vector3 p{point_cloud(i, 0), point_cloud(i, 1), point_cloud(i, 2)};
        point_set.push_back(p);
    }
    return get_convex_hull(point_set);
}


std::tuple<ManifoldSurfaceMesh*, VertexPositionGeometry*>
get_convex_hull_mesh(Eigen::MatrixX3d point_cloud){
    size_t n = point_cloud.rows();
    size_t dim = 3;
    std::vector<Vector3> point_set;
    for (size_t i = 0; i < n; i++){
        Vector3 p{point_cloud(i, 0), point_cloud(i, 1), point_cloud(i, 2)};
        point_set.push_back(p);
    }
    return get_convex_hull_mesh(point_set);
}


std::tuple<std::vector<std::vector<size_t>>, std::vector<size_t>, std::vector<Vector3>>
get_convex_hull(std::vector<Vector3> point_set){
    const size_t num_points = point_set.size();
    const size_t dim = 3;
    char flags[64];
    // sprintf(flags, "qhull Qt");
    std::vector<double> pset_flat_vec(num_points*dim);
    for (size_t i = 0; i < num_points ; i++){
        Vector3 p = point_set[i];
        pset_flat_vec[3*i] = p.x;
        pset_flat_vec[3*i + 1] = p.y;
        pset_flat_vec[3*i + 2] = p.z;
    }
    coordT* data = new coordT[num_points * dim];
    std::copy(pset_flat_vec.begin(), pset_flat_vec.end(), data);
    Qhull qhull("n", 3, num_points, data, "Qt"); // was using QJ Q3 and // Qt Q3 // then Qt ; Tv for debug?
    
    // make data structures
    size_t my_count = 0, max_ind = 0;
    std::vector<size_t> indices;
    for (QhullVertex qvertex: qhull.vertexList()){
        my_count++;
        indices.push_back(qvertex.id()-1);
        if (qvertex.id() - 1 > max_ind)
            max_ind = qvertex.id() - 1;
    }
    if (qhull.vertexCount() > max_ind + 1)
        printf(" vert count is %d, my count: %d, max_ind: %d \n", qhull.vertexCount(), my_count, max_ind);
    std::vector<std::vector<size_t>> hull_faces;
    std::vector<Vector3> poses(qhull.vertexCount());  // 
    std::sort(indices.begin(), indices.end());
    std::map<size_t, size_t> re_indexing;
    for (size_t i = 0; i < indices.size(); i++){
        re_indexing[indices[i]] = i;
    }

    // tracing old vertices
    std::vector<size_t> old_indices(qhull.vertexCount());
    for (QhullVertex qvertex: qhull.vertexList()){
        size_t old_ind = qvertex.point().id();
        old_indices[re_indexing[qvertex.id() - 1]] = old_ind;
        // printf("old index %d new index- %d\n", old_ind, re_indexing[qvertex.id() - 1]);
    }
    for (QhullFacet qfacet: qhull.facetList().toStdVector()){
        std::vector<size_t> hull_face;
        // std::cout << qfacet << "\n";
        for (QhullVertex qvertex: qfacet.vertices().toStdVector()){
            // size_t id = qvertex.id() - 1; // qhull is 1-based for some reason
            size_t id = re_indexing[qvertex.id() - 1]; // qhull is 1-based for some reason
            hull_face.push_back(id);
            
            // position
            Vector3 pos;
            const realT *c= qvertex.point().coordinates();
            for(int k= qvertex.dimension(); k--; ){
                // std::cout << " "<< k << " " << ; // QH11010 FIX: %5.2g
                if (k == 2) pos.x = *c++;
                if (k == 1) pos.y = *c++;
                if (k == 0) pos.z = *c++;
            }
            poses[id] = pos;
        }
        hull_faces.push_back(hull_face);
    }
    
    // build and return mesh
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;
    SurfaceMesh* surf_mesh = new SurfaceMesh(hull_faces);
    surf_mesh->greedilyOrientFaces();

    // make sure of outward orientation; can't just make qhull do it ??!
    // since mesh is already oriented and is convex, only need to check one face
    Face f0 = surf_mesh->face(0);
    Vertex v0 = f0.halfedge().tailVertex(),
           v1 = f0.halfedge().tipVertex(),
           v2 = f0.halfedge().next().tipVertex();
    Vector3 p0 = poses[v0.getIndex()],
            p1 = poses[v1.getIndex()],
            p2 = poses[v2.getIndex()];
    Vector3 p_other = poses[f0.halfedge().twin().next().tipVertex().getIndex()];
    if (dot(p_other - p0, cross(p1 - p0, p2 - p1)) > 0.){ 
        surf_mesh->invertOrientation(f0);
        surf_mesh->greedilyOrientFaces(); // will start with f0!
    }
    surf_mesh->compress();

    return {surf_mesh->getFaceVertexList(), old_indices, poses};
}


std::tuple<std::vector<std::vector<size_t>>, std::vector<size_t>, std::vector<Vector3>>
get_convex_hull(VertexData<Vector3> point_set){
    const size_t num_points = point_set.size();
    std::vector<Vector3> point_set_vec;
    for (size_t i = 0; i < num_points ; i++){
        Vector3 p = point_set[i];
        point_set_vec.push_back(p);
    }
    return get_convex_hull(point_set_vec);
}



std::tuple<ManifoldSurfaceMesh*, VertexPositionGeometry*>
get_convex_hull_mesh(VertexData<Vector3> point_set) {
    const size_t num_points = point_set.size();
    std::vector<Vector3> point_set_vec;
    for (size_t i = 0; i < num_points ; i++){
        Vector3 p = point_set[i];
        point_set_vec.push_back(p);
    }
    return get_convex_hull_mesh(point_set_vec);
}


std::tuple<ManifoldSurfaceMesh*, VertexPositionGeometry*>
get_convex_hull_mesh(std::vector<Vector3> point_set) {
    std::vector<std::vector<size_t>> hull_faces; 
    std::vector<size_t> hull_vertex_mapping;
    std::vector<Vector3> hull_poses; // redundant, but helps with keeping this function clean
    std::tie(hull_faces, hull_vertex_mapping, hull_poses) = get_convex_hull(point_set);
    ManifoldSurfaceMesh* hull_mesh = new ManifoldSurfaceMesh(hull_faces);
    VertexPositionGeometry* hull_geometry = new VertexPositionGeometry(*hull_mesh);
    for (Vertex hull_v: hull_mesh->vertices())
        hull_geometry->inputVertexPositions[hull_v] = hull_poses[hull_v.getIndex()];
    return std::tuple<ManifoldSurfaceMesh*, VertexPositionGeometry*>(hull_mesh, hull_geometry);
}



std::tuple<ManifoldSurfaceMesh*, VertexPositionGeometry*> 
get_mesh_for_convex_set(Eigen::MatrixX3d convex_point_cloud){
    std::vector<std::vector<size_t>> hull_faces; 
    std::vector<size_t> hull_vertex_mapping;
    std::vector<Vector3> hull_poses; // redundant, but helps with keeping this function clean
    std::tie(hull_faces, hull_vertex_mapping, hull_poses) = get_convex_hull(convex_point_cloud);
    
    // re-indexing
    std::vector<std::vector<size_t>> org_index_hull_faces; 
    for (std::vector<size_t> face: hull_faces){
        std::vector<size_t> org_face;
        for (size_t v: face){
            org_face.push_back(hull_vertex_mapping[v]);
        }
        org_index_hull_faces.push_back(org_face);
    }
    ManifoldSurfaceMesh *hull_mesh = new ManifoldSurfaceMesh(org_index_hull_faces);
    VertexPositionGeometry* hull_geometry = new VertexPositionGeometry(*hull_mesh);
    for (Vertex hull_v: hull_mesh->vertices())
        hull_geometry->inputVertexPositions[hull_v] = vec_to_GC_vec3(convex_point_cloud.row(hull_v.getIndex()));
    return std::tuple<ManifoldSurfaceMesh*, VertexPositionGeometry*>(hull_mesh, hull_geometry);
}



Vector3 project_back_into_hull(VertexPositionGeometry *hull_geometry, Vector3 p){
    // find closest face on hull
    
    Face closest_outward_face;
    double min_face_dist  = std::numeric_limits<double>::infinity();
    bool is_outside_flag = false;
    for (Face f: hull_geometry->mesh.faces()){ // check if face-projectable; while keeping track of closest vertex 
        Halfedge curr_he  = f.halfedge(),
                 first_he = f.halfedge();
        // assume outward normals; which is inward w.r.t. p
        Vector3 f_normal = hull_geometry->faceNormal(f);
        double face_dist = dot(f_normal, p - hull_geometry->inputVertexPositions[f.halfedge().vertex()]);
        if (face_dist >= 0) // outside the hull
            is_outside_flag = true;
        if (face_dist < 0) // inside the hull; at the wrong side of the face
            continue;
        bool face_is_projectable = true;
        while (true){ // checking if face-projectable
            Vertex v1 = curr_he.tailVertex(), v2 = curr_he.tipVertex(), v3 = curr_he.next().tipVertex();
            Vector3 A = hull_geometry->inputVertexPositions[v1], 
                    B = hull_geometry->inputVertexPositions[v2], 
                    C = hull_geometry->inputVertexPositions[v3];
            // assume outward normals; which is inward w.r.t. p
            Vector3 N_PAB = cross(A - B, p - A);
            if (dot(f_normal, N_PAB) >= 0) { // not face-projectable on this face
                face_is_projectable = false;
                break; // go to next face
            }
            curr_he = curr_he.next();
            if (curr_he == first_he)
                break;
        }
        if (face_is_projectable && face_dist < min_face_dist){
            min_face_dist = face_dist;
            closest_outward_face = f;
        }
    }

    if (!is_outside_flag)
        return p;
    
    // is outside the hull; use closest face if projectable
    if (min_face_dist != std::numeric_limits<double>::infinity()){ // use closest face
        // printf("face projected\n");
        return p - min_face_dist * hull_geometry->faceNormal(closest_outward_face) * 1.01;
    }

    // is outside the hull; not face-projectable; use closest edge
    double min_edge_dist = std::numeric_limits<double>::infinity();
    Edge closest_edge;
    Vector3 closest_edge_ortho_p;
    
    for (Edge e: hull_geometry->mesh.edges()){ // check if vertex-projectable
        Vector3 A = hull_geometry->inputVertexPositions[e.firstVertex()], 
                B = hull_geometry->inputVertexPositions[e.secondVertex()];
        Vector3 PB = B - p,
                PA = A - p,
                AB = B - A;
        bool edge_is_projectable = dot(PB, AB) >= 0 && dot(PA, -AB) >= 0; // similar to singular edge condition in rolling dynamics
        if (!edge_is_projectable)
            continue;
        Vector3 ortho_p = PB - AB*dot(AB, PB)/dot(AB,AB);
        double edge_dist =  ortho_p.norm();
        if (edge_dist < min_edge_dist){
            closest_edge = e;
            min_edge_dist = edge_dist;
            closest_edge_ortho_p = ortho_p;
        }
    }
    if (min_edge_dist != std::numeric_limits<double>::infinity()){ // no face-projectable face ; use closest projectable edge
        // printf("edge projected\n");
        return p + closest_edge_ortho_p * 1.01;
    }


    // not edge-projectable; use closest vertex
    double min_vertex_dist = std::numeric_limits<double>::infinity();
    Vertex closest_vertex;
    for (Vertex v: hull_geometry->mesh.vertices()){ 
        Vector3 V = hull_geometry->inputVertexPositions[v];
        double vertex_dist = (V - p).norm();
        if (vertex_dist < min_vertex_dist){
            closest_vertex = v;
            min_vertex_dist = vertex_dist;
        }
    }

    printf("vertex projected\n");
    return p + (hull_geometry->inputVertexPositions[closest_vertex] - p) * 1.01;
    
}