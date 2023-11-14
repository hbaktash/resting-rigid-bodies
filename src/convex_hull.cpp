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


// std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>> 
std::tuple<std::vector<std::vector<size_t>>, std::vector<size_t>, std::vector<Vector3>>
get_convex_hull(VertexData<Vector3> point_set){
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
    Qhull qhull("n", 3, num_points, data, "QJ Q3");
    
    // make data structures
    size_t my_count = 0, max_ind = 0;
    std::vector<size_t> indices;
    for (QhullVertex qvertex: qhull.vertexList()){
        my_count++;
        indices.push_back(qvertex.id()-1);
        if (qvertex.id() - 1 > max_ind)
            max_ind = qvertex.id() - 1;
    }
    if (qhull.vertexCount() > max_ind)
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

    // make sure of outward orientation; can't qhull just do it ?????!
    Face f0 = surf_mesh->face(0);
    Vertex v0 = f0.halfedge().tailVertex(),
           v1 = f0.halfedge().tipVertex(),
           v2 = f0.halfedge().next().tipVertex();
    Vector3 p0 = poses[v0.getIndex()],
            p1 = poses[v1.getIndex()],
            p2 = poses[v2.getIndex()];
    Vector3 p_other = poses[f0.halfedge().twin().next().tipVertex().getIndex()]; // damn you qhull!!!!
    if (dot(p_other - p0, cross(p1 - p0, p2 - p1)) > 0.){
        surf_mesh->invertOrientation(f0);
        surf_mesh->greedilyOrientFaces(); // will start with f0!
    }
    surf_mesh->compress();

    return {surf_mesh->getFaceVertexList(), old_indices, poses};
    // mesh = std::unique_ptr<ManifoldSurfaceMesh>(new ManifoldSurfaceMesh(surf_mesh->getFaceVertexList())); //  surf_mesh->toManifoldMesh()
    // // mesh.reset(new ManifoldSurfaceMesh(hull_faces));
    // geometry = std::unique_ptr<VertexPositionGeometry>(new VertexPositionGeometry(*mesh));
    // for (Vertex v : mesh->vertices()) {
    //     // Use the low-level indexers here since we're constructing
    //     (*geometry).inputVertexPositions[v] = poses[v.getIndex()];
    // }
    // return std::make_tuple(std::move(mesh), std::move(geometry));
}