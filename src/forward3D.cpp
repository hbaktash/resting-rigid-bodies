#include "forward3D.h"

Forward3DSolver::Forward3DSolver(ManifoldSurfaceMesh* inputMesh_, VertexPositionGeometry* inputGeo_,
                                 Vector3 inputG_, bool concave_input){
    inputMesh = inputMesh_;
    inputGeometry = inputGeo_;
    if (concave_input){
        std::vector<std::vector<size_t>> hull_faces; 
        std::vector<size_t> hull_vertex_mapping;
        std::vector<Vector3> hull_poses; // redundant, but helps with keeping this function clean
        std::tie(hull_faces, hull_vertex_mapping, hull_poses) = get_convex_hull(inputGeo_->inputVertexPositions);
        hullMesh = new ManifoldSurfaceMesh(hull_faces);
        hullGeometry = new VertexPositionGeometry(*hullMesh);
        org_hull_indices = VertexData<size_t>(*hullMesh);
        for (Vertex hull_v: hullMesh->vertices()){
            hullGeometry->inputVertexPositions[hull_v] = hull_poses[hull_v.getIndex()];
            org_hull_indices[hull_v] = hull_vertex_mapping[hull_v.getIndex()];
        }
    }
    else {
        hullMesh = inputMesh_;
        hullGeometry = inputGeo_;
        org_hull_indices = VertexData<size_t>(*hullMesh);
        for (Vertex hull_v: hullMesh->vertices())
            org_hull_indices[hull_v] = hull_v.getIndex();
    }
    // hullMesh = inputMesh_;
    // hullGeometry = inputGeo_;
    G = inputG_;
}


Forward3DSolver::Forward3DSolver(Eigen::MatrixX3d point_cloud, Eigen::Vector3d _G, bool is_convex){
    G = vec2vec3(_G);

    std::vector<std::vector<size_t>> hull_faces; 
    std::vector<size_t> hull_vertex_mapping;
    std::vector<Vector3> hull_poses; // redundant, but helps with keeping this function clean
    std::tie(hull_faces, hull_vertex_mapping, hull_poses) = get_convex_hull(point_cloud);
    
    if (!is_convex){
        hullMesh = new ManifoldSurfaceMesh(hull_faces);
        hullGeometry = new VertexPositionGeometry(*hullMesh);
        org_hull_indices = VertexData<size_t>(*hullMesh);
        for (Vertex hull_v: hullMesh->vertices()){
            hullGeometry->inputVertexPositions[hull_v] = hull_poses[hull_v.getIndex()];
            org_hull_indices[hull_v] = hull_vertex_mapping[hull_v.getIndex()];
        }
        inputMesh = hullMesh;
        inputGeometry = hullGeometry;
    }
    else {
        if (hull_poses.size() != point_cloud.rows()){ // make sure input was actually convex
            std::cout << "hull_poses size: " << hull_poses.size() << std::endl;
            std::cout << "point_cloud size: " << point_cloud.rows() << std::endl;
            throw std::logic_error(" Input set was not convex; hull poses size mismatch");
        }
        // re-indexing
        std::vector<std::vector<size_t>> org_index_hull_faces; 
        for (std::vector<size_t> face: hull_faces){
            std::vector<size_t> org_face;
            for (size_t v: face){
                org_face.push_back(hull_vertex_mapping[v]);
            }
            org_index_hull_faces.push_back(org_face);
        }
        hullMesh = new ManifoldSurfaceMesh(org_index_hull_faces);
        hullGeometry = new VertexPositionGeometry(*hullMesh);
        org_hull_indices = VertexData<size_t>(*hullMesh);
        for (Vertex hull_v: hullMesh->vertices()){
            hullGeometry->inputVertexPositions[hull_v] = vec_to_GC_vec3(point_cloud.row(hull_v.getIndex()));
            org_hull_indices[hull_v] = hull_v.getIndex(); // trivial assignment
        }
        inputMesh = hullMesh;
        inputGeometry = hullGeometry;
    }
    
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


// stabilizable := the normal from G can touch the ground
// stable := the normal from G falls withing the element
bool Forward3DSolver::vertex_is_stablizable(Vertex v){
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
    Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
    Vector3 p1 = hullGeometry->inputVertexPositions[v1], 
            p2 = hullGeometry->inputVertexPositions[v2];
    if (dot(G-p1, p2-p1) <= 0) return v1;  // <Gv1v2 is open 
    else if (dot(G-p2, p1-p2) <= 0) return v2;  // <Gv2v1 is open
    else return Vertex(); // INVALID ind for singular edges!
}

bool Forward3DSolver::edge_is_stable(Edge e){
    return next_rolling_vertex(e).getIndex() == INVALID_IND;
}



bool Forward3DSolver::edge_is_stablizable(Edge e){
    // finding the orthogonal vector from G onto AB 
    Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
    Vector3 A = hullGeometry->inputVertexPositions[v1], B = hullGeometry->inputVertexPositions[v2];
    Vector3 ortho_g =  point_to_segment_normal(G, A, B); //GB - AB*dot(AB, GB)/dot(AB,AB);

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
    // iterate over all edges of the face and check the plane angle made with G
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
        Face next_face;
        // printf("found angles\n");
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


void Forward3DSolver::face_to_next(Face f){
    if (face_is_stable(f)){
        // std::cout << "   -/-/-/-  at a stable face " << f.getIndex() << ": " << f.halfedge().tailVertex().getIndex() << ", " << f.halfedge().tipVertex().getIndex() << ", " << f.halfedge().next().tipVertex().getIndex() << "  -\\-\\-\\- \n";
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
        
        if (rolling_through_e0 && e0_is_singular) {
            Face next_face = curr_he.twin().face();
            curr_f = next_face;
            curr_v = Vertex();
            curr_e = Edge();
            curr_g_vec = hullGeometry->faceNormal(next_face);
            break;
        }
        if ((rolling_through_e0 || rolling_through_e1) && !e0_is_singular && !e1_is_singular) {
            // std::cout << "got to vertex!!\n";
            vertex_to_next(v1);
            break;
        }   
        // go to next he
        curr_he = curr_he.next();
        if (curr_he == first_he){
            throw std::logic_error("Face should be either stable or something should happen in the previous loop!\n");
        }
    }
}


void Forward3DSolver::next_state(bool verbose){
    if (stable_state){
        return;
    }
    if (curr_v.getIndex() != INVALID_IND){
        if (verbose){
            std::cout << " STATUS: at vertex " << curr_v.getIndex() << "\n";
            std::cout << "         curr gvec " << curr_g_vec << "\n";
        }
        Vertex old_v = curr_v;
        vertex_to_next(curr_v);
        if (curr_v == old_v){
            std::cerr << "Vertex " << curr_v.getIndex() << " is stable?!!??!\n";
            stable_state = true;
            return;
        }
    }
    else if (curr_e.getIndex() != INVALID_IND){
        if (verbose){
            std::cout << " STATUS: at edge " << curr_e.getIndex() << ": " << curr_e.firstVertex().getIndex() << ", " << curr_e.secondVertex().getIndex() << "\n";
        }
        Edge old_e = curr_e;
        edge_to_next(curr_e);
        if (curr_e == old_e){
            std::cerr << "Edge " << curr_e.getIndex() << " is stable?\n";
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
        std::cerr << " $$$ initialize the contact first! $$$" << std::endl;
    }
}


std::vector<Vector3> Forward3DSolver::quasi_static_drop(Vector3 initial_orientation){
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
    bool verbose = false;
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
    }
}

void Forward3DSolver::build_face_last_faces(){
    ///
    face_last_face = FaceData<Face>(*hullMesh, Face());
    for (Face f: hullMesh->faces()){
        if (face_last_face[f].getIndex() != INVALID_IND)
         continue;
        Face temp_f = f;
        std::vector<Face> faces_history;
        faces_history.push_back(temp_f);
        while (face_next_face[temp_f] != temp_f){
            temp_f = face_next_face[temp_f];
            faces_history.push_back(temp_f);
        }
        for (Face stored_f: faces_history) {
            face_last_face[stored_f] = temp_f;
        }
    }
}

// just call all the pre-compute initializations; not face-last-face
void Forward3DSolver::initialize_pre_computes(){
    compute_vertex_stabilizablity();
    compute_edge_stable_normals();
    build_face_next_faces();
    build_face_last_faces();
    inputGeometry->refreshQuantities();
    hullGeometry->refreshQuantities();
}


void Forward3DSolver::print_precomputes(){
    printf( " --------- face last faces %d --------- \n", hullMesh->nFaces());
    for (Face f: hullMesh->faces()){
        printf("face %d -> %d\n", f.getIndex(), face_last_face[f].getIndex());
    }
    printf( " --------- edge next vertex --------- \n");
    for (Edge e: hullMesh->edges()){
        printf("edge %d -> %d\n", e.getIndex(), edge_next_vertex[e].getIndex());
    }
    printf( " --------- vertex stability --------- \n");
    for (Vertex v: hullMesh->vertices()){
        if (vertex_is_stablizable(v))
            printf("vertex %d -> %d\n", v.getIndex(), vertex_is_stabilizable[v]);
    }
}

