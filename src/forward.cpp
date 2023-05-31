#include "forward.h"


ForwardSolver::ForwardSolver(ManifoldSurfaceMesh* inputMesh_, VertexPositionGeometry* inputGeo_,
                             Vector3 inputG_){
    mesh = inputMesh_;
    geometry = inputGeo_;
    G = inputG_;
}


void ForwardSolver::center_geometry(){
    Vector3 center{0.,0.,0.};
    for(Vertex v: mesh->vertices()){
        center += geometry->inputVertexPositions[v];
    }
    center /= (double) mesh->nVertices(); 
    for(Vertex v: mesh->vertices()){
        geometry->inputVertexPositions[v] -= center;
        hullGeometry->inputVertexPositions -= center;
    }
}

void ForwardSolver::align_geometry(Vector3 dir){
    center_geometry();
    double theta = atan2(dir.y, dir.x);
    double to_rotate = 1.5*PI - theta;
    for (Vertex v: hullMesh->vertices()){
        Vector3 pos = hullGeometry->inputVertexPositions[v];
        Vector2 tmp_v2{pos.x, pos.y};
        tmp_v2 = tmp_v2.rotate(to_rotate);
        hullGeometry->inputVertexPositions[v] = {tmp_v2.x, tmp_v2.y, 0.};
    }
}

void ForwardSolver::find_contact(Vector3 initial_ori){
    // just rotating the g_vec
    Vertex contact_point = mesh->vertex(0);
    double max_inner_product = 0.; // smth should be positive, if we choose a ref inside the convex hull
    for (Vertex v: hullMesh->vertices()){
        Vector3 pos = hullGeometry->inputVertexPositions[v],
                Gv = pos - G; // could use any point inside convexHull for this (using G just cause we have it).
        if (dot(Gv, initial_ori) >= max_inner_product){
            contact_point = v;
            max_inner_product = dot(Gv, initial_ori);
        }
    }
    Vertex v = contact_point;
    // update state & initialize tracing
    curr_state.first = v;
    curr_state.second = v;
    current_g_vec = initial_ori;
}

Edge ForwardSolver::other_edge(Edge curr_e, Vertex tip_v){
    for (Edge adj_e: tip_v.adjacentEdges()){
        if (adj_e != curr_e && 
            adj_e.isBoundary()){
                return adj_e;
        }
    }
}

void ForwardSolver::build_next_edge_tracer(){
    next_falling_edge = EdgeData<Edge>(*hullMesh, Edge()); // Valid for boundary edges.
    for (Edge e: hullMesh->edges()){
        if (e.isBoundary()){
            Vertex v1 = e.firstVertex(),
                v2 = e.secondVertex();
            Vector3 p1 = hullGeometry->inputVertexPositions[v1], 
                    p2 = hullGeometry->inputVertexPositions[v2];
            Vector3 p1g = G - p1,
                    p2g = G - p2;
            if (dot(p1g, p2-p1) < 0) // will roll to v1's side
                next_falling_edge[e] = other_edge(e, v1);
            else if(dot(p2g, p1-p2) < 0) // will roll to v2's side
                next_falling_edge[e] = other_edge(e, v2);
            else // stable edge
                next_falling_edge[e] = e;
        }
    }

    // test 
    for (Edge e: hullMesh->edges()){
        printf(" Edge %d next edge: %d \n", e.getIndex(), next_falling_edge[e].getIndex());
    }
}

void ForwardSolver::compute_vertex_probabilities(){
    vertex_probabilities = VertexData<double>(*hullMesh, 0.);
    for (Vertex v: hullMesh->vertices()){
        Edge e1, e2;
        Vertex v1, v2;
        bool first = true;
        for (Edge e: v.adjacentEdges()){
            if (e.isBoundary() && first){
                e1 = e;
                first = false;
                v1 = e1.otherVertex(v);
            }
            else if (e.isBoundary()){
                e2 = e;
                v2 = e2.otherVertex(v);
            }

        }
        Vector3 vv1 = hullGeometry->inputVertexPositions[v1] - hullGeometry->inputVertexPositions[v],
                vv2 = hullGeometry->inputVertexPositions[v2] - hullGeometry->inputVertexPositions[v];
        double angle = acos(dot(vv1, vv2)/(norm(vv1)*norm(vv2)));
        std::cout<< "angle at "<< v.getIndex() << " is "<< angle*180/PI << " \n";
        vertex_probabilities[v] = 0.5 * (1. - angle/PI);
    }
}

void ForwardSolver::compute_initial_edge_probabilities(){
    initial_edge_probabilities = EdgeData<double>(*hullMesh, 0.);
    for (Edge e: hullMesh->edges()){
        if (e.isBoundary()){
            Vector3 p1 = hullGeometry->inputVertexPositions[e.firstVertex()],
                    p2 = hullGeometry->inputVertexPositions[e.secondVertex()];
            Vector3 gp1 = p1 - G,
                    gp2 = p2 - G;
            double angle = acos(dot(gp1, gp2)/(norm(gp1)*norm(gp2)));
            initial_edge_probabilities[e] = angle/(2*PI);
            // printf(" edge prob %f \n", angle/(2*PI));
        }
    }
}

void ForwardSolver::compute_final_edge_probabilities(){
    // initials already should be computed
    final_edge_probabilities = EdgeData<double>(initial_edge_probabilities); // EdgeData<double>(*hullMesh, 0.);
    for (Edge e: hullMesh->edges()){
        if (e.isBoundary()){
            Edge final_edge = e; // finding the root/sink of the DAG
            while (next_falling_edge[final_edge] != final_edge) {
                final_edge = next_falling_edge[final_edge];
            }
            // pour the probability from e to final_edge
            double current_prob = final_edge_probabilities[e];
            final_edge_probabilities[final_edge] += current_prob;
            final_edge_probabilities[e] -= current_prob;
        }
    }
    for (Edge e: hullMesh->edges()){
        if (e.isBoundary()){
            printf("Edge %d final prob: %f \n", e.getIndex(), final_edge_probabilities[e]);
        }
    }
}

void ForwardSolver::next_state(){
    if (curr_state.first == curr_state.second){
        // we are on a vertex
        // mesh.
    }
    else {
        // we are on an edge

    }
}
