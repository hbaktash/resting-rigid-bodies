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
    // for (Edge e: hullMesh->edges()){
    //     printf(" Edge %d next edge: %d \n", e.getIndex(), next_falling_edge[e].getIndex());
    // }
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
    printf("---- Final edge probs ---- \n");
    for (Edge e: hullMesh->edges()){
        if (e.isBoundary()){
            printf("Edge %d final prob: %f \n", e.getIndex(), final_edge_probabilities[e]);
        }
    }
}


void ForwardSolver::find_contact(Vector3 initial_ori){ // 0-based vector, need to add G_offset
    // just rotating the g_vec
    Vertex contact_point = mesh->vertex(0);
    double max_inner_product = 0.; // smth should be positive, if we choose a ref inside the convex hull
    for (Vertex v: hullMesh->vertices()){
        Vector3 pos = hullGeometry->inputVertexPositions[v];
        Vector3 Gv = pos - G; // could use any point inside convexHull for this (using G just cause we have it).
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

    // also build next edge tracer for next() calls
}


void ForwardSolver::next_state(){
    std::pair<Vertex, Vertex> next_state;
    Vector3 next_g_vec;
    if (curr_state.first == curr_state.second){
        // we are on a vertex
        Vertex v = curr_state.first;
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
        Vector3 Gp  = hullGeometry->inputVertexPositions[v] - G , 
                Gp1 = hullGeometry->inputVertexPositions[v1]- G ,
                Gp2 = hullGeometry->inputVertexPositions[v2]- G ; // these 3 vectors are sufficient to determine next falling element
        Vector3 cross_p = cross(current_g_vec, Gp),
                cross_p1 = cross(current_g_vec, Gp1),
                cross_p2 = cross(current_g_vec, Gp2);
        next_state.first = v; // dah
        // v1, v2 on different sides
        if (dot(cross_p, cross_p1) >= 0 && dot(cross_p, cross_p2) <= 0) // fall on to v-v2
            next_state.second = v2;
        else if (dot(cross_p, cross_p1) <= 0 && dot(cross_p, cross_p2) >= 0) // fall on to v-v1
            next_state.second = v1;
        else { // both v1,v2 in one side
            Vector3 cross_p1p = cross(Gp1, Gp),
                    cross_p2p = cross(Gp2, Gp);
            if (dot(cross_p1p, cross_p) >= 0 && dot(cross_p2p, cross_p) <= 0)
                next_state.second = v1;
            else if (dot(cross_p2p, cross_p) >= 0 && dot(cross_p1p, cross_p) <= 0)
                next_state.second = v2;
            else 
                printf("what the hell happend?\n");   
        }
        assert(next_state.first != next_state.second);
        Vector3 edge_vec = hullGeometry->inputVertexPositions[next_state.second] - hullGeometry->inputVertexPositions[next_state.first],
                z = {0., 0., 1.};
        next_g_vec = cross(edge_vec, z);
        if (dot(next_g_vec, Gp) >= 0) // wrong orientaion
            next_g_vec = -next_g_vec;
        
        // update things
        curr_state = next_state;
        current_g_vec = next_g_vec;
    }
    else {
        // we are on an edge
        Edge current_edge = hullMesh->connectingEdge(curr_state.first, curr_state.second);
        Edge next_edge = next_falling_edge[current_edge];
        // update things
        // update state (skipping rolling on vertices)
        curr_state.first = next_edge.firstVertex();
        curr_state.second = next_edge.secondVertex();
        // update g-vec
        Vertex v = curr_state.first;
        Vector3 Gp = G - hullGeometry->inputVertexPositions[v];
        Vector3 edge_vec = hullGeometry->inputVertexPositions[curr_state.second] - hullGeometry->inputVertexPositions[curr_state.first],
                z = {0., 0., 1.};
        next_g_vec = cross(edge_vec, z);
        if (dot(next_g_vec, Gp) >= 0) // wrong orientaion
            next_g_vec = -next_g_vec;
        current_g_vec = next_g_vec;
    }
}


Edge ForwardSolver::simulate_toss(Vector3 initial_ori){
    find_contact(initial_ori);
    next_state();
    Edge current_edge = hullMesh->connectingEdge(curr_state.first, curr_state.second);
    while (next_falling_edge[current_edge] != current_edge) {
        next_state();
        current_edge = hullMesh->connectingEdge(curr_state.first, curr_state.second);
    }
    return current_edge;
    // return Edge();
}


void ForwardSolver::empirically_build_probabilities(int sample_count){
    empirical_final_probabilities = EdgeData<double>(*hullMesh, 0.);
    int succ_sample_count = 0;
    for (int i = 0; i < sample_count; i++){
        double x = randomReal(-1, 1),
               y = randomReal(-1, 1);
        double norm = x*x + y*y;
        if (norm <= 1){ // rejecting norm>=1 samples
            succ_sample_count += 1;
            Vector3 initial_ori = {x/sqrt(norm), y/sqrt(norm), 0.};
            Edge final_edge = simulate_toss(initial_ori);
            empirical_final_probabilities[final_edge] += 1;
        }
    }
    for (Edge e: hullMesh->edges()){

        empirical_final_probabilities[e] /= (double)succ_sample_count;
        if (empirical_final_probabilities[e] != 0)
            printf("empirical probability Edge %d: %f\n", e.getIndex(), empirical_final_probabilities[e]);
    }
}
