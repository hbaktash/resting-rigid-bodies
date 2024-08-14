
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
#include "MS_complex.h"

QuasiDynamicSolver::QuasiDynamicSolver(ManifoldSurfaceMesh* inputMesh_, VertexPositionGeometry* inputGeo_,
                             Vector3 inputG_, bool concave_input){
    // inputMesh = inputMesh_;
    // inputGeometry = inputGeo_;
    // if (concave_input){
    //     update_convex_hull();
    // }
    // else {
    //     hullMesh = inputMesh_;
    //     hullGeometry = inputGeo_;
    //     trivial_initialize_index_trackers(); //
    // }
    // // hullMesh = inputMesh_;
    // // hullGeometry = inputGeo_;
    // G = inputG_;
    // updated = false;
}


QuasiDynamicSolver::QuasiDynamicSolver(Eigen::MatrixX3d _point_cloud, Eigen::Vector3d _G){
    std::vector<std::vector<size_t>> hull_faces; 
    std::vector<size_t> hull_vertex_mapping;
    std::vector<Vector3> hull_poses; // redundant, but helps with keeping this function clean
    auto [_inputMesh, _inputGeometry] = get_convex_hull_mesh(point_cloud);
    input_points = _point_cloud.vec;
    hullMesh = inputMesh;
    hullGeometry = inputGeometry;
    G = Vector3{_G[0], _G[1], _G[2]};
    updated = false;
}

// initialize state holders
void Forward3DSolver::initialize_state(Vertex curr_v_, Edge curr_e_, Face curr_f_, Vector3 curr_g_vec_){
    curr_v = curr_v_;
    curr_e = curr_e_;
    curr_f = curr_f_;
    curr_g_vec = curr_g_vec_;
    stable_state = face_is_stable(curr_f); // null-ness is checked inside the function
}


//
void Forward3DSolver::set_G(Vector3 new_G){
    this->G = new_G;
    updated = false;
}

void Forward3DSolver::set_uniform_G(){
    std::pair<Vector3, double> G_V_pair = find_center_of_mass(*inputMesh, *inputGeometry);
    set_G(G_V_pair.first);
    volume = G_V_pair.second;
}
//
Vector3 Forward3DSolver::get_G(){
    return G; 
}

void Forward3DSolver::trivial_initialize_index_trackers(){
    // hull_indices = .resize(hullMesh->nVertices());
    // interior_indices = Vector<size_t>::Zero(0);
    org_hull_indices = VertexData<size_t>(*hullMesh);
    on_hull_index = VertexData<size_t>(*inputMesh, INVALID_IND);    
    for (Vertex hull_v: hullMesh->vertices())
        org_hull_indices[hull_v] = hull_v.getIndex();
    for (Vertex org_v: inputMesh->vertices())
        on_hull_index[org_v] = org_v.getIndex();
}

void Forward3DSolver::update_hull_index_arrays(){
    hull_indices.resize(hullMesh->nVertices());
    interior_indices.resize(inputMesh->nVertices() - hullMesh->nVertices());
    for (Vertex v: hullMesh->vertices()){
        hull_indices[v.getIndex()] = org_hull_indices[v];
    }
    size_t cnt = 0;
    for (Vertex v: inputMesh->vertices()){
        if (on_hull_index[v] == INVALID_IND){
            // printf("updating interior v %d \n", v.getIndex());
            interior_indices[cnt++] = v.getIndex();
        }
    }
}

void Forward3DSolver::update_convex_hull(bool with_projection){
    // printf(" -- updating convex hull -- \n");
    if (!with_projection || first_hull){ // just update and take the new hull
        std::vector<std::vector<size_t>> hull_faces; 
        std::vector<size_t> hull_vertex_mapping;
        std::vector<Vector3> hull_poses; // redundant, but helps with keeping this function clean
        std::tie(hull_faces, hull_vertex_mapping, hull_poses) = get_convex_hull(inputGeometry->inputVertexPositions);

        hullMesh = new ManifoldSurfaceMesh(hull_faces);
        hullGeometry = new VertexPositionGeometry(*hullMesh);
        org_hull_indices = VertexData<size_t>(*hullMesh);
        on_hull_index = VertexData<size_t>(*inputMesh, INVALID_IND);
        hull_indices.resize(hullMesh->nVertices());
        
        for (Vertex v: hullMesh->vertices()){
            Vector3 new_pos = hull_poses[v.getIndex()];
            hullGeometry->inputVertexPositions[v] = new_pos;
            inputGeometry->inputVertexPositions[hull_vertex_mapping[v.getIndex()]] = new_pos;
            org_hull_indices[v] = hull_vertex_mapping[v.getIndex()];
            on_hull_index[hull_vertex_mapping[v.getIndex()]] = v.getIndex();
            hull_indices[v.getIndex()] = org_hull_indices[v];
        }

        first_hull = false;
    }
    else { // hull gets smaller always in vertex count
        // update hull pos first
        for (Vertex v: hullMesh->vertices()){
            hullGeometry->inputVertexPositions[v] = inputGeometry->inputVertexPositions[org_hull_indices[v]];
        }
        // getting the hull of the current (possibly updated) hull
        std::vector<std::vector<size_t>> hull_hull_faces; 
        std::vector<size_t> hull_hull_vertex_mapping;
        std::vector<Vector3> hull_hull_poses; // redundant, but helps with keeping this function clean
        std::tie(hull_hull_faces, hull_hull_vertex_mapping, hull_hull_poses) = get_convex_hull(hullGeometry->inputVertexPositions); 

        // update index trackers
        hullMesh = new ManifoldSurfaceMesh(hull_hull_faces);
        hullGeometry = new VertexPositionGeometry(*hullMesh);
        
        VertexData<size_t> new_org_hull_indices = VertexData<size_t>(*hullMesh);
        on_hull_index = VertexData<size_t>(*inputMesh, INVALID_IND); // no need to preserve the previous mapping
        hull_indices.resize(hullMesh->nVertices());
        for (Vertex hull_v: hullMesh->vertices()){
            // index updates
            new_org_hull_indices[hull_v] = org_hull_indices[hull_hull_vertex_mapping[hull_v.getIndex()]]; // m1(m2(x))
            on_hull_index[new_org_hull_indices[hull_v]] = hull_v.getIndex();
            hull_indices[hull_v.getIndex()] = new_org_hull_indices[hull_v];
            // geo updates
            Vector3 new_pos = hull_hull_poses[hull_v.getIndex()];
            hullGeometry->inputVertexPositions[hull_v] = new_pos;
            // inputGeometry->inputVertexPositions[new_org_hull_indices[v]] = new_pos;
        }
        org_hull_indices = new_org_hull_indices;

        polyscope::registerSurfaceMesh("pre-projection mesh", inputGeometry->inputVertexPositions,
                                                              inputMesh->getFaceVertexList())->setEnabled(false);
        // projection step
        for (Vertex v: inputMesh->vertices()){
            if (on_hull_index[v] == INVALID_IND){ // interior vertex
                Vector3 new_pos = project_back_into_hull(hullGeometry, inputGeometry->inputVertexPositions[v]);
                if (new_pos != inputGeometry->inputVertexPositions[v]){
                    // printf("projected vertex %d\n", v.getIndex());
                }
                inputGeometry->inputVertexPositions[v] = new_pos;
            }
        }
    }
    // printf(" -- convex hull updated -- \n");
    
    // updating index vectors
    update_hull_index_arrays();
}

// redo the indices
void Forward3DSolver::update_hull_points_correspondence(VertexData<Vector3> new_hull_points, VertexData<Vector3> old_points){
    VertexData<bool> assigned(*inputMesh, false),
                     best_match_taken(*hullMesh, false);
    on_hull_index = VertexData<size_t>(*inputMesh, INVALID_IND);
    org_hull_indices = VertexData<size_t>(*hullMesh, INVALID_IND);
    for (Vertex v: hullMesh->vertices()){ // hull points
        Vector3 hull_p = new_hull_points[v];
        double min_dist = 1e10;
        Vertex best_v = Vertex();
        for (Vertex other_v: inputMesh->vertices()){
            double tmp_dist = (hull_p - old_points[other_v]).norm();
            if (tmp_dist <= min_dist){
                if (!assigned[other_v]){
                    min_dist = tmp_dist;
                    best_v = other_v;
                    best_match_taken[v] = false;
                }
                else 
                    best_match_taken[v] = true;
            }
        }
        assigned[best_v] = true;
        
        assert(best_v.getIndex() != INVALID_IND);
        // Vertex prev_org_vertex = inputMesh->vertex(org_hull_indices[v]);
        // printf("hull v %d assigned to %d\n", v.getIndex(), best_v.getIndex());
        // if (best_match_taken[v])
        //     printf(" :((( best match taken\n");
        // on_hull_index[prev_org_vertex] = INVALID_IND;
        on_hull_index[best_v] = v.getIndex();
        org_hull_indices[v] = best_v.getIndex();
    }
    update_hull_index_arrays();
}

// height from G to ground after contact
double Forward3DSolver::height_function(Vector3 ground_normal){
    Vertex contact_point = hullMesh->vertex(0);
    double max_inner_product = 0.; // smth should be positive, if we choose a ref inside the convex hull
    for (Vertex v: hullMesh->vertices()){
        Vector3 pos = hullGeometry->inputVertexPositions[v];
        Vector3 Gv = pos - G; // could use any point inside the convexHull (using G just just cause we have it).
        if (dot(Gv, ground_normal) >= max_inner_product){
            contact_point = v;
            max_inner_product = dot(Gv, ground_normal);
        }
    }
    return max_inner_product;
}


// find initial contact and basically refresh the current state
void Forward3DSolver::find_contact(Vector3 initial_ori){
    // just find the max inner product with the g_vec 
    Vertex contact_point = hullMesh->vertex(0);
    double max_inner_product = 0.; // smth should be positive, if we choose a ref inside the convex hull
    for (Vertex v: hullMesh->vertices()){
        Vector3 pos = hullGeometry->inputVertexPositions[v];
        Vector3 Gv = pos - G; // could use any point inside the convexHull (using G just just cause we have it).
        if (dot(Gv, initial_ori) >= max_inner_product){
            contact_point = v;
            max_inner_product = dot(Gv, initial_ori);
        }
    }

    // update state & initialize tracing
    curr_v = contact_point;
    curr_e = Edge();
    curr_f = Face();

    curr_g_vec = initial_ori;

    stable_state = false;

    // for vector field visuals
    Vector3 p = hullGeometry->inputVertexPositions[contact_point];
    Vector3 G_proj = project_on_plane(G, p, curr_g_vec.normalize()); // project G onto the ground plane going through P
    initial_roll_dir = (G_proj - p);
}


// just the Gaussian curvature
void Forward3DSolver::compute_vertex_probabilities(){
    vertex_probabilities = VertexData<double>(*hullMesh, 0.);
    for (Vertex v: hullMesh->vertices()){
        vertex_probabilities[v] = hullGeometry->vertexGaussianCurvature(v)/(4*PI);
    }
}

// stabilizable := the normal from G can touch the ground
// stable := the normal from G falls withing the element
bool Forward3DSolver::vertex_is_stablizable(Vertex v){
    if (updated)
        return vertex_is_stabilizable[v];
    Vector3 Gp = hullGeometry->inputVertexPositions[v] - G;
    double gp_norm = dot(Gp, Gp);
    for (Vertex other_v: hullMesh->vertices()){
        if (v != other_v){
            Vector3 Gp2 = hullGeometry->inputVertexPositions[other_v] - G;
            if (dot(Gp2, Gp) > gp_norm)
                return false;
        }
    }
    return true;
}


Vertex Forward3DSolver::next_rolling_vertex(Edge e){
    if (updated)
        return edge_next_vertex[e];
    Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
    Vector3 p1 = hullGeometry->inputVertexPositions[v1], 
            p2 = hullGeometry->inputVertexPositions[v2];
    if (dot(G-p1, p2-p1) <= 0) return v1;  // <Gv1v2 is open 
    else if (dot(G-p2, p1-p2) <= 0) return v2;  // <Gv2v1 is open
    else return Vertex(); // INVALID ind for singular edges!
}

bool Forward3DSolver::edge_is_stable(Edge e){
    if (updated)
        return edge_next_vertex[e].getIndex() == INVALID_IND;
    return next_rolling_vertex(e).getIndex() == INVALID_IND;
}



bool Forward3DSolver::edge_is_stablizable(Edge e){
    // finding the orthogonal vector from G onto AB 
    Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
    Vector3 A = hullGeometry->inputVertexPositions[v1], B = hullGeometry->inputVertexPositions[v2];
    Vector3 ortho_g =  point_to_segment_normal(G, A, B); //GB - AB*dot(AB, GB)/dot(AB,AB);
    
    // debugger
    tmp_test_vec = ortho_g;
    // checking if other vertices are more aligned with orhto_G
    double ortho_g_norm = dot(ortho_g, ortho_g);
    for (Vertex other_v: hullMesh->vertices()){
        if (other_v != v1 && other_v != v2){
            Vector3 GP = hullGeometry->inputVertexPositions[other_v] - G;
            if (dot(GP,ortho_g) > ortho_g_norm){
                return false;
            }
        }
    }
    return true;
}


bool Forward3DSolver::face_is_stable(Face f){
    if (f.getIndex() == INVALID_IND)
        return false;
    if (updated)
        return face_next_face[f] == f;
    // need to check all the wideness of dihedral angles made on each edge
    // iterate over all edges of the face and check the dihedral angle
    Halfedge curr_he = f.halfedge(),
             first_he = f.halfedge();
    while (true){
        Vertex v1 = curr_he.tailVertex(), v2 = curr_he.tipVertex(), v3 = curr_he.next().tipVertex();
        Vector3 A = hullGeometry->inputVertexPositions[v1], 
                B = hullGeometry->inputVertexPositions[v2], 
                C = hullGeometry->inputVertexPositions[v3];
        // assume outward normals!
        Vector3 N_ABC = hullGeometry->faceNormal(f),// cross(B - A, C - B),
                N_GAB  = cross(A - B, G - A);
        if (dot(N_ABC, N_GAB) > 0)
            return false;
        curr_he = curr_he.next();
        if (curr_he == first_he)
            break;
    }
    return true;
}


void Forward3DSolver::vertex_to_next(Vertex v){
    // printf(" doing vertex to next\n");
    Vector3 p = hullGeometry->inputVertexPositions[v],
            unit_g_vec = curr_g_vec.normalize();
    Vector3 Gp = p - G;
    Vector3 G_proj = project_on_plane(G, p, unit_g_vec); // project G onto the ground plane going through P
    Vector3 rotation_plane_normal = cross(-Gp, G_proj - p).normalize(); // this plane doesn't change until smth hits the ground
    
    Vector3 pG_proj = G_proj - p; // the projected PG vector

    // find the vertex with the least angle made with projected PG
    double max_dot_cos = -1.; // since convex, max angle will be PI
    Vertex best_v = Vertex(); // next vertex that hits the ground
    // NOTE: from here projections refer to projection on the rotation axis plane
    Vector3 best_p2_proj; // to help with finding next g_vec after the loop
    for (Vertex v2: v.adjacentVertices()){
        Vector3 p2 = hullGeometry->inputVertexPositions[v2];
        Vector3 p2_proj = project_on_plane(p2, p, rotation_plane_normal); // project neigh vertices onto the rotation axis plane
        Vector3 pp2_proj = p2_proj - p;
        double tmp_dot_cos = dot(pG_proj, pp2_proj)/(norm(pG_proj)*norm(pp2_proj));
        if (tmp_dot_cos > max_dot_cos){ // robust?
            max_dot_cos = tmp_dot_cos;
            best_v = v2;
            best_p2_proj = p2_proj;
        }
    }
    
    // update current state
    curr_e = hullMesh->connectingEdge(v, best_v);
    curr_v = Vertex();
    curr_f = Face();
    // computing next g_vec
    Vector3 unit_next_pp = (p - best_p2_proj).normalize();
    Vector3 next_g_vec = Gp - unit_next_pp * dot(unit_next_pp, Gp); 
    curr_g_vec = next_g_vec.normalize();
}


void Forward3DSolver::edge_to_next(Edge e){
    Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
    Vector3 p1 = hullGeometry->inputVertexPositions[v1], p2 = hullGeometry->inputVertexPositions[v2];
    // if (e.getIndex() == 177){
    //     printf("doing edge to next\n");
    //     std::cout << " p1: " << p1 << "\n p2:" << p2 << "\n";
    //     std::cout << "   G: "<<G << "\n";
    //     std::cout << "   G + currG: "<<G + curr_g_vec << "\n";
    // }
    if (dot(G - p1, p2-p1) < 0.){ // if the Gp1p2 angle is wide
        vertex_to_next(v1); // rolls to the v1 vertex
    }
    else if (dot(G - p2, p1-p2) < 0.){ // if the Gp2p1 angle is wide
        vertex_to_next(v2); // rolls to the v2 vertex
    }
    else { // rolls to a neighboring face
        // printf("rolling to face from edge\n");
        Vector3 Gp1 = p1 - G,
                unit_g_vec = curr_g_vec.normalize();
        Vector3 G_proj = G + unit_g_vec * dot(Gp1, unit_g_vec), // project G on the ground
                rotation_plane_normal = (p2 - p1).normalize();
        // future projections will be onto the rotation plane

        // find immediate neighbors of v1 on the two neigh faces; 
        // assuming planar polygonal faces
        Halfedge he = e.halfedge();
        Vertex va = he.prevOrbitFace().tailVertex(), // va is on the he.face()  
               vb = he.twin().next().tipVertex();    // vb is on the twin.face()
        
        // project va, vb, projected_G onto the rotation plane
        // the plane passes through "p1" with the "rotation axis" as the normal
        Vector3 pa = hullGeometry->inputVertexPositions[va],
                pb = hullGeometry->inputVertexPositions[vb];
        Vector3 G_proj_proj = project_on_plane(G_proj, p1, rotation_plane_normal),
                pa_proj = project_on_plane(pa, p1, rotation_plane_normal),
                pb_proj = project_on_plane(pb, p1, rotation_plane_normal);
        double pa_angle_pre_acos = dot(pa_proj - p1, G_proj_proj - p1)/(norm(pa_proj - p1)*norm(G_proj_proj - p1)),
               pb_angle_pre_acos = dot(pb_proj - p1, G_proj_proj - p1)/(norm(pb_proj - p1)*norm(G_proj_proj - p1));
        double pa_angle = acos(pa_angle_pre_acos),
               pb_angle = acos(pb_angle_pre_acos);
        // if (e.getIndex() == 177){
        //     std::cout << "   DEBUG dot values: " << dot(pa_proj - p1, G_proj_proj - p1) << " " << dot(pb_proj - p1, G_proj_proj - p1) << "\n";
        //     std::cout << "   norms a         : " << norm(pa_proj - p1) << " " << norm(G_proj_proj - p1) << "\n";
        //     std::cout << "   norms b         : " << norm(pb_proj - p1) << " " << norm(G_proj_proj - p1) << "\n";
        //     std::cout << "   DEBUG pre acos v: " << pa_angle_pre_acos << " " << pb_angle_pre_acos << "\n";
        //     std::cout << "   acos debug: " << pa_angle << " " << pb_angle << "\n";
        // }
        Face next_face;
        // printf("found angles\n");
        // if (pa_angle <= pb_angle){ // face containing va is next
        if (pa_angle_pre_acos >= pb_angle_pre_acos){
            curr_f = he.face();
            curr_v = Vertex();
            curr_e = Edge();
        }
        else { // face containing vb is next
            curr_f = he.twin().face();
            curr_v = Vertex();
            curr_e = Edge();
        }
        // always assuming outward normals
        curr_g_vec = hullGeometry->faceNormal(curr_f);
    }
}

// assuming face is not stable
void Forward3DSolver::face_to_next(Face f){
    if (face_is_stable(f)){
        // printf("   -/-/-/-  at a stable face %d: %d, %d, %d  -\\-\\-\\- \n",f.getIndex(), f.halfedge().tailVertex().getIndex(), f.halfedge().tipVertex().getIndex(), f.halfedge().next().tipVertex().getIndex());
        stable_state = true;
        return;
    }
    Vector3 unit_g_vec = curr_g_vec.normalize();
    if (dot(G - hullGeometry->inputVertexPositions[f.halfedge().vertex()], curr_g_vec) >= 0)
        printf("____G is outside_____!! \n");
    // project G_proj on the same plane as the face plane
    assert(hullGeometry->faceNormal(f) == curr_g_vec); // assume outward normals
    Vector3 G_proj = project_on_plane(G, hullGeometry->inputVertexPositions[f.halfedge().tailVertex()], unit_g_vec);
    Halfedge curr_he = f.halfedge(),
             first_he = f.halfedge();
    while (true) {
        // printf("at he %d\n", curr_he.getIndex());
        // will only check rolling to v0v1 at this iteration; v2 is only used for checking rolling onto v1
        Vertex v1 = curr_he.tipVertex(),
               v0 = curr_he.tailVertex(),
               v2 = curr_he.next().tipVertex();
        Vector3 p1 = hullGeometry->inputVertexPositions[v1],
                p0 = hullGeometry->inputVertexPositions[v0],
                p2 = hullGeometry->inputVertexPositions[v2];
        // first check if we are on the "exterior" side of v0v1
        
        bool rolling_through_e0 = false;
        Vector3 N_Gp0p1 = cross(p0 - G, p1 - p0);
        if (dot(unit_g_vec, N_Gp0p1) <= 0)
            rolling_through_e0 = true;

        bool rolling_through_e1 = false;
        Vector3 N_Gp1p2 = cross(p1 - G, p2 - p1);
        if (dot(unit_g_vec, N_Gp1p2) <= 0)
            rolling_through_e1 = true;
        
        bool e0_is_singular = edge_is_stable(curr_he.edge()),
             e1_is_singular = edge_is_stable(curr_he.next().edge());
        
        // // Old version
        // Vector3 on_plane_edge_normal = cross(unit_g_vec, p0 - p1);
        // // use p2 as a point on the interior side
        // double p2_dot = dot(p2-p0, on_plane_edge_normal),
        //        G_proj_dot = dot(G_proj - p0, on_plane_edge_normal);
        // double angle_Gv0v1 = acos(dot(G_proj - p0, p1 - p0)/(norm(G_proj - p0)*norm(p1 - p0))),
        //        angle_Gv1v0 = acos(dot(G_proj - p1, p0 - p1)/(norm(G_proj - p1)*norm(p0 - p1))),
        //        angle_Gv1v2 = acos(dot(G_proj - p1, p2 - p1)/(norm(G_proj - p1)*norm(p2 - p1)));
        // if (p2_dot * G_proj_dot <= 0 && angle_Gv0v1 <= PI/2. && angle_Gv1v0 <= PI/2.){ // G' is on the exterior side; and angles are acute
        if (rolling_through_e0 && e0_is_singular) {
            Face next_face = curr_he.twin().face();
            curr_f = next_face;
            curr_v = Vertex();
            curr_e = Edge();
            curr_g_vec = hullGeometry->faceNormal(next_face);
            // printf("got to face!!\n");
            break;
        }
        // if (angle_Gv1v0 >= PI/2. && angle_Gv1v2 >= PI/2.){
        if ((rolling_through_e0 || rolling_through_e1) && !e0_is_singular && !e1_is_singular) {
            // printf("got to vertex!!\n");
            vertex_to_next(v1);
            break;
        }   
        // go to next he
        curr_he = curr_he.next();
        // std::cout << "at face "  <<  f.getIndex() << " " << v0.getIndex() << " " << v1.getIndex() << " " << v2.getIndex() << "\n";
        // std::cout <<    "poses: " << p0 << "\n \t\t" << p1 << "\n \t\t" << p2 << "\n";
        // std::cout <<    "G loc: " << G << "\n";
        // std::cout << "quantities: " << rolling_through_e0 << " "<< rolling_through_e1 << " " << e0_is_singular << " " << e1_is_singular << "\n";
        if (curr_he == first_he){
            throw std::logic_error("Face should be either stable or something should happen in the previous loop!\n");
        }
    }
}


void Forward3DSolver::next_state(bool verbose){
    if (stable_state){
        // printf(" &&&& already at a stable state! &&&&\n");
        return;
    }
    int status_check = (curr_v.getIndex() != INVALID_IND) + (curr_e.getIndex() != INVALID_IND) + (curr_f.getIndex() != INVALID_IND);
    // printf("test for bool addition %d\n", status_check);
    assert(status_check <= 1); // either not initiated (0) or only one valid (1)
    if (curr_v.getIndex() != INVALID_IND){
        if (verbose){
            printf(" STATUS: at vertex %d\n", curr_v.getIndex());
            std::cout << "         curr gvec " << curr_g_vec << "\n";
        }
        Vertex old_v = curr_v;
        vertex_to_next(curr_v);
        if (curr_v == old_v){
            printf("Vertex %d is stable?!!??!\n", curr_v.getIndex());
            stable_state = true;
            return;
        }
    }
    else if (curr_e.getIndex() != INVALID_IND){
        if (verbose){
            printf(" STATUS: at edge %d: %d, %d\n", curr_e.getIndex(), curr_e.firstVertex().getIndex(), curr_e.secondVertex().getIndex());
            // std::cout << "         curr gvec " << curr_g_vec << "\n";
            // std::cout << " edge vertex poses: " << hullGeometry->inputVertexPositions[curr_e.firstVertex()] << " -- " << hullGeometry->inputVertexPositions[curr_e.secondVertex()] << "\n ";
            // std::cout << " center of mass pose: " << G << "\n";
        }
        Edge old_e = curr_e;
        edge_to_next(curr_e);
        if (curr_e == old_e){
            printf("Edge %d is stable?!!??!\n", curr_e.getIndex());
            stable_state = true;
            return;
        }
    }
    else if (curr_f.getIndex() != INVALID_IND){
        Face old_f = Face(curr_f); // TODO: do I have to do this?
        Halfedge he = curr_f.halfedge();
        if (verbose){
            printf(" STATUS: at face %d: %d, %d, %d, ..\n", old_f.getIndex(), he.tailVertex().getIndex(), he.tipVertex().getIndex(), he.next().tipVertex().getIndex());
            std::cout << "         curr gvec " << curr_g_vec << "\n";
        }
        face_to_next(curr_f);
        if (curr_f == old_f){
            stable_state = true;
            return;
        }
    }
    else {
        printf(" $$$ initialize the contact first! $$$\n");
    }
}


std::vector<Vector3> Forward3DSolver::snail_trail_log(Vector3 initial_orientation){
    std::vector<Vector3> trail;
    translation_log.empty();
    vertex_log.empty();
    edge_log.empty();
    face_log.empty();
    trail.push_back(initial_orientation);
    find_contact(initial_orientation);
    vertex_log.push_back(curr_v);
    edge_log.push_back(curr_e);
    face_log.push_back(curr_f);    
    while (!stable_state){
        vertex_log.push_back(curr_v);
        edge_log.push_back(curr_e);
        face_log.push_back(curr_f);
        next_state();
        trail.push_back(curr_g_vec);
        if (stable_state)
            break;
    }
    return trail;
}


Face Forward3DSolver::final_touching_face(Vector3 initial_ori){
    find_contact(initial_ori);
    // printf(" Found contact! \n");
    while(true){
        next_state();
        if (stable_state)
            break;
    }
    // assert(curr_f.getIndex() != INVALID_IND);
    return curr_f;
}



/// pre computes


// is 0 if the normal is unreachable; and if non-singular: the normal doesnt fall on the edge (edge too short)
void Forward3DSolver::compute_edge_stable_normals(){
    edge_next_vertex = EdgeData<Vertex>(*hullMesh, Vertex()); 
    Vector3 zero_vec({0.,0.,0.});
    edge_stable_normal = EdgeData<Vector3>(*hullMesh, zero_vec);
    for (Edge e: hullMesh->edges()){
        Vertex next_vertex = next_rolling_vertex(e);
        if (next_vertex.getIndex() == INVALID_IND){ // singular edge
            Vector3 A = hullGeometry->inputVertexPositions[e.firstVertex()],
                    B = hullGeometry->inputVertexPositions[e.secondVertex()];
            edge_stable_normal[e] = point_to_segment_normal(G, A, B).normalize();
        }
        edge_next_vertex[e] = next_vertex;
    }
}

// pre-compute vertex gaussian curvatures
void Forward3DSolver::compute_vertex_gaussian_curvatures(){
    vertex_gaussian_curvature = VertexData<double>(*hullMesh, 0.);
    for (Vertex v: hullMesh->vertices()){
        vertex_gaussian_curvature[v] = gaussian_curvature(v, *hullGeometry);
        // printf(" Gaussian Curvature at %d, is %f\n", v.getIndex(), vertex_gaussian_curvature[v]);
    }
}

// whether the stable point can be reached or not
void Forward3DSolver::compute_vertex_stabilizablity(){
    vertex_is_stabilizable = VertexData<bool>(*hullMesh, false);
    vertex_stable_normal = VertexData<Vector3>(*hullMesh);
    for (Vertex v: hullMesh->vertices()){
        if (vertex_is_stablizable(v))
            vertex_is_stabilizable[v] = true;
        vertex_stable_normal[v] = (hullGeometry->inputVertexPositions[v] - G).normalize();
    }
}


// deterministically find the next rolling face 
void Forward3DSolver::build_face_next_faces(){
    face_next_face = FaceData<Face>(*hullMesh);
    bool verbose = false; //f.getIndex() == 812 || f.getIndex() == 105; //
    // if (verbose)
    // printf("building face next faces\n");
    for (Face f: hullMesh->faces()){
        if (verbose) printf("at face %d\n", f.getIndex());
        initialize_state(Vertex(), Edge(), f, hullGeometry->faceNormal(f)); // assuming outward normals
        next_state(verbose); // could roll to an edge
        size_t count = 0;
        while (curr_f.getIndex() == INVALID_IND){
            count++;
            next_state(verbose);
        } // terminates when it gets to the next face
        if (verbose)
            printf(" last stable face was %d\n", curr_f.getIndex());
        face_next_face[f] = curr_f;
        // printf(" fnf %d -> %d\n", f.getIndex(), curr_f.getIndex());
    }
}

void Forward3DSolver::build_face_last_faces(){
    ///
    if (!updated)
        build_face_next_faces();
    face_last_face = FaceData<Face>(*hullMesh, Face());
    for (Face f: hullMesh->faces()){
        if (face_last_face[f].getIndex() != INVALID_IND)
         continue;
        // printf(" at face %d\n", f.getIndex());
        Face temp_f = f;
        std::vector<Face> faces_history;
        faces_history.push_back(temp_f);
        while (face_next_face[temp_f] != temp_f){
            temp_f = face_next_face[temp_f];
            // printf(" next face is %d\n", temp_f.getIndex());
            faces_history.push_back(temp_f);
        } // terminates when it gets to the next face
        for (Face stored_f: faces_history) {
            face_last_face[stored_f] = temp_f;
        }
        // printf(" ++ final face %d \n", temp_f.getIndex());
    }
}

// just call all the pre-compute initializations; not face-last-face
void Forward3DSolver::initialize_pre_computes(){
    // printf("precomputes:\n");
    // printf("  vertex stability:\n");
    compute_vertex_stabilizablity();
    // printf("  vertex gauss curvature:\n");
    compute_vertex_gaussian_curvatures();
    // printf("  edge stability:\n");
    compute_edge_stable_normals();
    // printf("  building face next faces:\n");
    build_face_next_faces(); // 
    // printf("  building face last faces:\n");
    build_face_last_faces();
    // printf("precomputes done!\n");
    inputGeometry->refreshQuantities();
    // printf("heree\n");
    hullGeometry->refreshQuantities();
    // printf("hereee\n");
    updated = true;
}
