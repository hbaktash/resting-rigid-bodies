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

#include "boundary_tools.h"


size_t BoundaryNormal::counter = 0;

// constructor
BoundaryNormal::BoundaryNormal(Vector3 _normal)
    : index(counter++){
        normal = _normal;
}

void BoundaryNormal::add_neighbor(BoundaryNormal* _neigh) {
    neighbors.push_back(_neigh);
}

// constructor
BoundaryBuilder::BoundaryBuilder(Forward3DSolver *forward_solver_){
    forward_solver = forward_solver_;
    mesh = forward_solver->hullMesh;
    geometry = forward_solver->hullGeometry;
}


 void BoundaryBuilder::build_boundary_normals(){
    vertex_boundary_normals = VertexData<BoundaryNormal*>(*mesh, nullptr);
    edge_boundary_normals = EdgeData<BoundaryNormal*>(*mesh, nullptr);

    // 
    std::list<Edge> bfs_list;
    EdgeData<bool> edge_visited(*mesh, false);
    VertexData<bool> vertex_visited(*mesh, false);
    // 
    forward_solver->build_face_last_faces(); // calls face next face within
    for (Edge e: mesh->edges()){
        if (forward_solver->edge_is_singular[e]){
            Face f1 = e.halfedge().face(),
                 f2 = e.halfedge().twin().face();
            if (forward_solver->face_last_face[f1] != forward_solver->face_last_face[f2]){ 
                // proved: at least one singular edge like this must exist
                // TODO: assert that stable normal falls inside the edge arc 
                Vector3 stable_edge_normal = forward_solver->edge_stable_normal[e];
                BoundaryNormal *new_boundary_normal = new BoundaryNormal(stable_edge_normal);
                edge_boundary_normals[e] = new_boundary_normal;

                bfs_list.push_back(e);
            }   
        }
    }
    // find them all..
    while (bfs_list.size() > 0) {
        Edge tmp_e = bfs_list.front();
        bfs_list.pop_front();
        edge_visited[tmp_e] = true;
        for (Vertex v: tmp_e.adjacentVertices()) {
            if (forward_solver->vertex_is_stabilizable[v]) {
                if(!vertex_visited[v]){ // normal is created if visited
                    BoundaryNormal *new_boundary_normal = new BoundaryNormal(forward_solver->vertex_stable_normal[v]);
                    vertex_boundary_normals[v] = new_boundary_normal;
                    vertex_visited[v] = true;
                }
                vertex_boundary_normals[v]->add_neighbor(edge_boundary_normals[tmp_e]);
                edge_boundary_normals[tmp_e]->add_neighbor(vertex_boundary_normals[v]);
            } // this vertex will be reached from the other side as well; no need to handle that here
            else { // unstable vertex, only one other boundary edge exists on this vertex patch
                if (vertex_visited[v]) continue; // skipping since vertex is unstable and already handled fully
                for (Edge e: v.adjacentEdges()){
                    if (!edge_visited[e]){
                        Face f1 = e.halfedge().face(),
                             f2 = e.halfedge().twin().face();
                        if (forward_solver->face_last_face[f1] != forward_solver->face_last_face[f2]){
                            Vector3 tmp_e_bnd_normal = edge_boundary_normals[tmp_e]->normal,
                                    vertex_stable_normal = forward_solver->vertex_stable_normal[v],
                                    f1_normal = geometry->faceNormal(f1),
                                    f2_normal = geometry->faceNormal(f2);
                            Vector3 e_bnd_normal = intersect_arc_ray_with_arc(vertex_stable_normal, tmp_e_bnd_normal, f1_normal, f2_normal);
                            assert(e_bnd_normal.norm() > 1e-7); // should intersect
                            // new boundary normal
                            BoundaryNormal *new_boundary_normal = new BoundaryNormal(e_bnd_normal);
                            edge_boundary_normals[e] = new_boundary_normal;
                            edge_boundary_normals[tmp_e]->add_neighbor(new_boundary_normal);
                            new_boundary_normal->add_neighbor(edge_boundary_normals[tmp_e]);
                            // for later iterations
                            bfs_list.push_back(e);
                            vertex_visited[v] = true;
                        }
                    }
                }
                assert(vertex_visited[v] == true);
            }
        }
    }
    
 }