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
        host_e = Edge();
        host_v = Vertex();
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
    edge_boundary_normals  = EdgeData<std::vector<BoundaryNormal*>>(*mesh);
    
    BoundaryNormal::counter = 0;
    // 
    std::vector<Edge> terminal_edges;
    // 
    printf("  buidling face-last-face\n");
    forward_solver->build_face_last_faces(); // calls face next face within
    printf("  finding terminal edges \n");
    for (Edge e: mesh->edges()){
        if (forward_solver->edge_next_vertex[e].getIndex() == INVALID_IND){
            Face f1 = e.halfedge().face(),
                 f2 = e.halfedge().twin().face();
            if (forward_solver->face_last_face[f1] != forward_solver->face_last_face[f2]){ 
                // proved: at least one singular edge like this must exist
                // TODO: assert that stable normal falls inside the edge arc 
                Vector3 stable_edge_normal = forward_solver->edge_stable_normal[e];
                BoundaryNormal *new_boundary_normal = new BoundaryNormal(stable_edge_normal);
                edge_boundary_normals[e].push_back(new_boundary_normal);
                
                new_boundary_normal->f1 = forward_solver->face_last_face[f1];
                new_boundary_normal->f2 = forward_solver->face_last_face[f2];
                new_boundary_normal->host_e = e;
                
                terminal_edges.push_back(e);
            }
        }
    }

    // for quick assignment of face-boundary-loops
    face_attraction_boundary = FaceData<std::vector<BoundaryNormal*>>(*mesh);
    face_region_area = FaceData<double>(*mesh, 0.);
    // back-flow them all..
    printf("  back-flowing terminal edges \n");
    for (Edge e: terminal_edges){
        // printf("- starting at terminal edge: %d\n", e.getIndex());
        assert(edge_boundary_normals[e].size() == 1);
        BoundaryNormal *bnd_normal = edge_boundary_normals[e].front();
        for (Vertex v: {e.firstVertex(), e.secondVertex()}){
            Vector3 tmp_normal = bnd_normal->normal,
                    f1_normal = geometry->faceNormal(bnd_normal->f1),
                    f2_normal = geometry->faceNormal(bnd_normal->f2),
                    v_normal = forward_solver->vertex_stable_normal[v];
            // Vector3 imm_f1_normal = geometry->faceNormal(e.halfedge().face()), // immediate face neighbors
            //         imm_f2_normal = geometry->faceNormal(e.halfedge().twin().face());
            
            double f1_area_sign = dot(f1_normal, cross(v_normal, tmp_normal)) >= 0 ? 1. : -1.;
            if (forward_solver->vertex_is_stabilizable[v])
                flow_back_boundary_on_edge(bnd_normal, Edge(), v, 
                                               f1_area_sign, f1_normal, f2_normal);
            else {
                for (Edge neigh_e: v.adjacentEdges()){
                    if (neigh_e != e){
                        // face_attraction_boundary[bnd_normal->f1].push_back(bnd_normal);
                        // face_attraction_boundary[bnd_normal->f2].push_back(bnd_normal);
                        flow_back_boundary_on_edge(bnd_normal, neigh_e, v, 
                                                f1_area_sign, f1_normal, f2_normal);
                    }
                }
            }
        }
    }
    for (Face f: mesh->faces()){
        if (face_region_area[f] > 0)
            face_region_area[f] /= (4.*PI); // sum to one
    }}




// recursively follow the boundary curve to a source
void BoundaryBuilder::flow_back_boundary_on_edge(BoundaryNormal* bnd_normal, Edge src_e, Vertex common_vertex,
                                                 double f1_area_sign, Vector3 f1_normal, Vector3 f2_normal){
    FaceData<std::vector<std::tuple<BoundaryNormal*, BoundaryNormal*, double>>> face_chain_area;
    // bnd_normal has to be in boundary normals of dest_e; won't assert tho for better performance
    Vertex v = common_vertex; // given as argument for better performance
    // printf("back-flowing to vertex %d\n", v.getIndex());
    Vector3 tmp_normal({0.,0.,0.}), next_normal({0.,0.,0.});
    if (forward_solver->vertex_is_stabilizable[v]){ // we are at a source
        BoundaryNormal *curr_vertex_boundary_normal = vertex_boundary_normal[v];
        if (vertex_boundary_normal[v] == nullptr){ // first time arriving at this stable vertex; create the boundary normal
            Vector3 stable_vertex_normal = forward_solver->vertex_stable_normal[v];
            curr_vertex_boundary_normal = new BoundaryNormal(stable_vertex_normal);
            vertex_boundary_normal[v] = curr_vertex_boundary_normal;
            curr_vertex_boundary_normal->host_v = v;
        }
        curr_vertex_boundary_normal->add_neighbor(bnd_normal);
        bnd_normal->add_neighbor(curr_vertex_boundary_normal);
        // handling face region boundary
        // face_attraction_boundary[bnd_normal->f1].push_back(curr_vertex_boundary_normal);
        // face_attraction_boundary[bnd_normal->f2].push_back(curr_vertex_boundary_normal);

        tmp_normal = bnd_normal->normal;
        next_normal = curr_vertex_boundary_normal->normal;
    }
    else {
        // vertex is not a source
        if (forward_solver->edge_next_vertex[src_e] == v){ // src_e is a source for this vertex
            Face f1 = src_e.halfedge().face(),
                f2 = src_e.halfedge().twin().face();
                
            // ** actually the following condition does not necessarily have to hold
            //    commented:
            // if (forward_solver->face_last_face[f1] == forward_solver->face_last_face[f2])
            //     return; // this source edge is not a source for this boundary normal 

            Vector3 vertex_stable_normal = forward_solver->vertex_stable_normal[v],
                    tmp_f1_normal = geometry->faceNormal(f1),
                    tmp_f2_normal = geometry->faceNormal(f2);
            Vector3 e_bnd_normal = intersect_arc_ray_with_arc(vertex_stable_normal, bnd_normal->normal, 
                                                            tmp_f1_normal, tmp_f2_normal);
            // printf("at tmp_e %d, %d. side v %d. e: %d, %d \n", dest_e.firstVertex().getIndex(), dest_e.secondVertex().getIndex(), v.getIndex(), e.firstVertex().getIndex(), e.secondVertex().getIndex());
            if(e_bnd_normal.norm() == 0.) 
                return; // not a source for the given bnd_normal
            // found the source normal


            BoundaryNormal *new_boundary_normal = new BoundaryNormal(e_bnd_normal);
            edge_boundary_normals[src_e].push_back(new_boundary_normal);
            bnd_normal->add_neighbor(new_boundary_normal);
            new_boundary_normal->add_neighbor(bnd_normal);
            new_boundary_normal->host_e = src_e;

            // region labels
            new_boundary_normal->f1 = bnd_normal->f1;
            new_boundary_normal->f2 = bnd_normal->f2;
            face_attraction_boundary[f1].push_back(new_boundary_normal);
            face_attraction_boundary[f2].push_back(new_boundary_normal);

            tmp_normal  = bnd_normal->normal,
            next_normal = new_boundary_normal->normal;
            
            // go with the back-flow
            Vertex next_v = src_e.otherVertex(v);
            if (forward_solver->vertex_is_stabilizable[next_v])
                flow_back_boundary_on_edge(new_boundary_normal, Edge(), next_v,
                                                f1_area_sign, f1_normal, f2_normal);
            else{
                for (Edge next_src_e: next_v.adjacentEdges()){
                    if (next_src_e != src_e){
                        // non_singularity and divisive-ness will be checked inside the function
                        flow_back_boundary_on_edge(new_boundary_normal, next_src_e, next_v,
                                                    f1_area_sign, f1_normal, f2_normal);
                    }
                }
            }
        }
    }
    // Vector3 dir_checker = ;
    if (tmp_normal.norm() == 0.)
        return;
    double curr_f1_alignment = dot(f1_normal, cross(next_normal, tmp_normal)) >= 0 ? 1. : -1.;
    double curr_f2_alignment = dot(f2_normal, cross(next_normal, tmp_normal)) >= 0 ? 1. : -1.;
    double f1_sign_change = f1_area_sign == curr_f1_alignment ? 1. : -1;
    double f2_sign_change = (-f1_area_sign == curr_f2_alignment) ? 1.: -1;
    // printf(" at v: %d f1, f2 sign change %f, %f \n", v.getIndex(),f1_sign_change, f2_sign_change);
    face_region_area[bnd_normal->f1] += 
                f1_sign_change * 
                    triangle_patch_area_on_sphere(f1_normal, bnd_normal->normal, next_normal);
    face_region_area[bnd_normal->f2] += 
                f2_sign_change * 
                    triangle_patch_area_on_sphere(f2_normal, bnd_normal->normal, next_normal);
    // TODO: take care of when f1,f2 on the same side when starting from saddle 
}

void BoundaryBuilder::print_area_of_boundary_loops(){
    printf(" Face probs:\n");
    for (Face f: mesh->faces()){
        if (face_region_area[f] > 0){
            // face_region_area[f] /= (4.*PI);
            printf(" f %d: %f\n", f.getIndex(), face_region_area[f]);
        }
    }
    // std::cout << face_region_area.toVector().transpose()<< "\n";
}