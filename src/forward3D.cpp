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
#include "forward3D.h"


// forward solver functions
Vector3 project_on_plane(Vector3 point, Vector3 offset, Vector3 normal){
    Vector3 unit_normal = normal.normalize();
    return point + unit_normal * dot(offset - point, unit_normal);
}

Vector3 point_to_segment_normal(Vector3 P, Vector3 A, Vector3 B){
    Vector3 PB = B - P,
            AB = B - A;
    Vector3 ortho_p = PB - AB*dot(AB, PB)/dot(AB,AB);
    return ortho_p;
}

Forward3DSolver::Forward3DSolver(ManifoldSurfaceMesh* inputMesh_, VertexPositionGeometry* inputGeo_,
                             Vector3 inputG_){
    hullMesh = inputMesh_;
    hullGeometry = inputGeo_;
    G = inputG_;
}

// initialize state holders
void Forward3DSolver::initialize_state(Vertex curr_v_, Edge curr_e_, Face curr_f_, Vector3 curr_g_vec_){
    curr_v = curr_v_;
    curr_e = curr_e_;
    curr_f = curr_f_;
    curr_g_vec = curr_g_vec_;
    stable_state = face_is_stable(curr_f); // null-ness is checked inside the function
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
    // need to check all the dihedral angles made on each edge
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
        if (dot(N_ABC, N_GAB) >= 0)
            return false;
        curr_he = curr_he.next();
        if (curr_he == first_he)
            break;
    }
    return true;
}


void Forward3DSolver::vertex_to_next(Vertex v){
    Vector3 p = hullGeometry->inputVertexPositions[v],
            unit_g_vec = curr_g_vec.normalize();
    Vector3 Gp = p - G;
    Vector3 G_proj = project_on_plane(G, p, unit_g_vec); // project G onto the ground plane going through P
    Vector3 rotation_plane_normal = cross(-Gp, G_proj - p).normalize(); // this plane doesn't change until smth hits the ground
    
    Vector3 pG_proj = G_proj - p; // the projected PG vector

    // find the vertex with the least angle made with projected PG
    double min_angle = PI; // since convex, max angle will be PI
    Vertex best_v = Vertex(); // next vertex that hits the ground
    // NOTE: from here projections refer to projection on the rotation axis plane
    Vector3 best_p2_proj; // to help with finding next g_vec after the loop
    for (Vertex v2: v.adjacentVertices()){
        Vector3 p2 = hullGeometry->inputVertexPositions[v2];
        Vector3 p2_proj = project_on_plane(p2, p, rotation_plane_normal); // project neigh vertices onto the rotation axis plane
        Vector3 pp2_proj = p2_proj - p;
        double tmp_angle = acos(dot(pG_proj, pp2_proj)/(norm(pG_proj)*norm(pp2_proj)));
        if (tmp_angle <= min_angle){
            min_angle = tmp_angle;
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
    curr_g_vec = next_g_vec;
}


void Forward3DSolver::edge_to_next(Edge e){
    Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
    Vector3 p1 = hullGeometry->inputVertexPositions[v1], p2 = hullGeometry->inputVertexPositions[v2];
    if (dot(G - p1, p2-p1) <= 0){ // if the Gp1p2 angle is wide
        vertex_to_next(v1); // rolls to the v1 vertex
    }
    else if (dot(G - p2, p1-p2) <= 0){ // if the Gp2p1 angle is wide
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
        double pa_angle = acos(dot(pa_proj - p1, G_proj_proj - p1)/(norm(pa_proj - p1)*norm(G_proj_proj - p1))),
               pb_angle = acos(dot(pb_proj - p1, G_proj_proj - p1)/(norm(pb_proj - p1)*norm(G_proj_proj - p1)));
        Face next_face;
        // printf("found angles\n");
        if (pa_angle <= pb_angle){ // face containing va is next
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
        Vector3 on_plane_edge_normal = cross(unit_g_vec, p0 - p1);
        // use p2 as a point on the interior side
        double p2_dot = dot(p2-p0, on_plane_edge_normal),
               G_proj_dot = dot(G_proj - p0, on_plane_edge_normal);
        double angle_Gv0v1 = acos(dot(G_proj - p0, p1 - p0)/(norm(G_proj - p0)*norm(p1 - p0))),
               angle_Gv1v0 = acos(dot(G_proj - p1, p0 - p1)/(norm(G_proj - p1)*norm(p0 - p1))),
               angle_Gv1v2 = acos(dot(G_proj - p1, p2 - p1)/(norm(G_proj - p1)*norm(p2 - p1)));
        if (p2_dot * G_proj_dot <= 0 && angle_Gv0v1 <= PI/2. && angle_Gv1v0 <= PI/2.){ // G' is on the exterior side; and angles are acute
            Face next_face = curr_he.twin().face();
            curr_f = next_face;
            curr_v = Vertex();
            curr_e = Edge();
            curr_g_vec = hullGeometry->faceNormal(next_face);
            break;
        }
        if (angle_Gv1v0 >= PI/2. && angle_Gv1v2 >= PI/2.){
            // printf("got to vertex!!\n");
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


void Forward3DSolver::next_state(){
    if (stable_state){
        // printf(" &&&& already at a stable state! &&&&\n");
        return;
    }
    int status_check = (curr_v.getIndex() != INVALID_IND) + (curr_e.getIndex() != INVALID_IND) + (curr_f.getIndex() != INVALID_IND);
    // printf("test for bool addition %d\n", status_check);
    assert(status_check <= 1); // either not initiated (0) or only one valid (1)
    if (curr_v.getIndex() != INVALID_IND){
        // printf(" STATUS: at vertex %d\n", curr_v.getIndex());
        Vertex old_v = curr_v;
        vertex_to_next(curr_v);
        if (curr_v == old_v){
            printf("Vertex %d is stable?!!??!\n", curr_v.getIndex());
            stable_state = true;
            return;
        }
    }
    else if (curr_e.getIndex() != INVALID_IND){
        // printf(" STATUS: at edge %d: %d, %d\n", curr_e.getIndex(), curr_e.firstVertex().getIndex(), curr_e.secondVertex().getIndex());
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
        // printf(" STATUS: at face %d: %d, %d, %d, ..\n", old_f.getIndex(), he.tailVertex().getIndex(), he.tipVertex().getIndex(), he.next().tipVertex().getIndex());
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
