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
    vertex_boundary_normal = VertexData<BoundaryNormal*>(*mesh, nullptr);
    edge_boundary_normals = EdgeData<std::vector<BoundaryNormal*>>(*mesh);

    // 
    std::vector<Edge> terminal_edges;
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
                edge_boundary_normals[e].push_back(new_boundary_normal);
                
                terminal_edges.push_back(e);
            }   
        }
    }
    // back-flow them all..
    for (Edge e: terminal_edges){
        assert(edge_boundary_normals[e].size() == 1);
        BoundaryNormal *bnd_normal = edge_boundary_normals[e].front();
        for (Vertex v: e.adjacentVertices()){
            if (forward_solver->vertex_is_stabilizable[v]){
                Vector3 stable_vertex_normal = forward_solver->vertex_stable_normal[v];
                BoundaryNormal *new_boundary_normal = new BoundaryNormal(stable_vertex_normal);
                new_boundary_normal->add_neighbor(bnd_normal);
                bnd_normal->add_neighbor(new_boundary_normal);
            }
            else {
                for (Edge neigh_e: v.adjacentEdges()){
                    if (neigh_e != e){
                        flow_back_boundary_on_edge(e, bnd_normal, neigh_e, v);
                    }
                }
            }
        }
    }
    
 }


void BoundaryBuilder::flow_back_boundary_on_edge(Edge dest_e, BoundaryNormal* bnd_normal, Edge src_e, Vertex common_vertex){
    // bnd_normal has to be in boundary normals of dest_e; won't assert tho for better performance
    
    if (forward_solver->edge_is_singular[src_e])
        return; // then src_e cannot really be a source 
    Vertex v = common_vertex; // given as argument for better performance
    Face f1 = src_e.halfedge().face(),
         f2 = src_e.halfedge().twin().face();
         
    if (forward_solver->face_last_face[f1] == forward_solver->face_last_face[f2])
        return; // this source edge is not a source for this boundary normal 

    Vector3 vertex_stable_normal = forward_solver->vertex_stable_normal[v],
            f1_normal = geometry->faceNormal(f1),
            f2_normal = geometry->faceNormal(f2);
    Vector3 e_bnd_normal = intersect_arc_ray_with_arc(vertex_stable_normal, bnd_normal->normal, 
                                                        f1_normal, f2_normal);
    // printf("at tmp_e %d, %d. side v %d. e: %d, %d \n", dest_e.firstVertex().getIndex(), dest_e.secondVertex().getIndex(), v.getIndex(), e.firstVertex().getIndex(), e.secondVertex().getIndex());
    if(e_bnd_normal.norm() == 0.) 
        return; // not a source for the given bnd_normal
    // found the source normal
    BoundaryNormal *new_boundary_normal = new BoundaryNormal(e_bnd_normal);
    edge_boundary_normals[src_e].push_back(new_boundary_normal);
    bnd_normal->add_neighbor(new_boundary_normal);
    new_boundary_normal->add_neighbor(bnd_normal);

    // follow the backflow
    Vertex next_v = src_e.otherVertex(v);
    for (Edge next_src_e: next_v.adjacentEdges()){
        if (next_src_e != src_e){
            // non_singularity and divisive-ness will be checked inside the function
            flow_back_boundary_on_edge(src_e, new_boundary_normal, next_src_e, next_v);
        }
    }
}