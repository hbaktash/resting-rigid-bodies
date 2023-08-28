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

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
// #include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
// #include "bullet3/examples/BasicExample.h"
#include "args/args.hxx"
#include "imgui.h"

#include "coloring.h"
#include "forward3D.h"
#include "mesh_factory.h"
#include "geometry_utils.h"
#include "visual_utils.h"
#include "markov_model.h"

// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Polyhedron_3.h>
// #include <CGAL/Surface_mesh.h>
// #include <CGAL/convex_hull_3.h>
// #include <vector>
// #include <fstream>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_ptr;
std::unique_ptr<VertexPositionGeometry> geometry_ptr;
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;
Vector3 G, // center of Mass
        initial_g_vec({0,-1,0}),
        default_face_color({0.99,0.99,0.99}),
        curr_face_color({0.1,0.87,0.1});

Forward3DSolver forwardSolver;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psInputMesh, *dummy_psMesh1, 
                       *dummy_psMesh2, *dummy_psMesh3, *dummy_psMesh_for_patches,
                       *dummy_forward_vis,
                       *coloredPsMesh;

polyscope::PointCloud *psG, // point cloud with single G
                      *curr_state_pt, *curr_g_vec_gm_pt,
                      *gauss_map_pc, *face_normals_pc,
                      *stable_face_normals_pc, *edge_equilibria_pc, *stabilizable_edge_pc, *stable_edge_pc,
                      *stable_vertices_gm_pc, *hidden_stable_vertices_gm_pc,
                      *raster_pc;

polyscope::SurfaceGraphQuantity* curr_state_segment;

float pt_cloud_radi_scale = 0.1,
      curve_radi_scale = 0.1,
      G_r = 0.,
      G_theta = 0.,
      G_phi = 0.,
      g_vec_theta = 0.,
      g_vec_phi = 0.;

float stable_edge_radi = 0.007,
      stablizable_edge_radi = 0.009,
      both_edge_radi = 0.013,
      pt_cloud_stablizable_radi = 0.03,
      face_normal_vertex_gm_radi = 0.03;
glm::vec3 stable_edge_color({0.2, 0.3, 0.3}),
          stabilizable_edge_color({0.05, 0.05, 0.05}),
          // stabilizable_edge_color({0.1, 0.4, 0.6}),
          both_edge_color({0.1, 0.1, 0.8});


// Gauss map stuff
double gm_distance = 2.1,
       gm_radi = 1.0;
Vector3 gm_shift({0., gm_distance, 0.}),
        colored_shift({gm_distance, gm_distance, 0.});
bool color_arcs = false,
     draw_unstable_edge_arcs = true,
     draw_stable_g_vec_for_unstable_edge_arcs = false,
     show_hidden_stable_vertex_normals = true,
     draw_stable_patches = false;

// arc stuff
float arc_curve_radi = 0.01;
glm::vec3 patch_arc_fancy_color({0.9,0.1,0.1});

int arcs_seg_count = 13,
    arc_counter = 0;

// raster image stuff
int sample_count = 8000; 
FaceData<Vector3> face_colors;
bool recolor_faces = true,
     real_time_raster = false,
     normalize_vecF = true;
glm::vec3 vecF_color({0.1, 0.1, 0.1});

// snail trail stuff
bool show_snail_trail = true;
polyscope::PointCloud *snail_trail_pc;
polyscope::SurfaceMesh *dummy_ps_mesh_for_snail_trail;
Vector3 old_g_vec, new_g_vec;
int snail_trail_dummy_counter = 0;

// Markov Chain stuff
RollingMarkovModel *markov_model;
polyscope::PointCloud *sudo_faces_pc;

// example choice
std::vector<std::string> all_polyhedra_items = {std::string("cube"), std::string("tet"), std::string("sliced tet"), std::string("Conway spiral 4")};
std::string all_polygons_current_item = "tet";
static const char* all_polygons_current_item_c_str = "tet";


void draw_stable_patches_on_gauss_map(){
  std::vector<std::vector<size_t>> dummy_face{{0,0,0}};
  dummy_psMesh_for_patches = polyscope::registerSurfaceMesh("dummy mesh for patch arcs", geometry->inputVertexPositions, dummy_face);
  if (draw_stable_patches){
    arc_counter = 0;
    for (Edge e: forwardSolver.hullMesh->edges()){
      for (Vertex v: e.adjacentVertices()){
        if (forwardSolver.edge_is_stable(e) && forwardSolver.edge_is_stablizable(e) &&
            forwardSolver.vertex_is_stablizable(v)){
          Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
          Vector3 A = forwardSolver.hullGeometry->inputVertexPositions[v1], B = forwardSolver.hullGeometry->inputVertexPositions[v2];
          Vector3 GB = B - G,
                  AB = B - A;
          Vector3 ortho_g = GB - AB*dot(AB, GB)/dot(AB,AB);
          Vector3 v_stable_vec = forwardSolver.hullGeometry->inputVertexPositions[v] - forwardSolver.G;
          draw_arc_on_sphere(normalize(ortho_g), normalize(v_stable_vec), gm_shift, gm_radi, arcs_seg_count, 100 + arc_counter, dummy_psMesh_for_patches, 1., 
                            color_arcs ? patch_arc_fancy_color : glm::vec3({-1.,0,0})); // fancy color if coloring edge arcs
          arc_counter++;
        }
      }
    }
  }
}

void draw_stable_vertices_on_gauss_map(){
  Vector3 shift = {0., gm_distance, 0.};
  std::vector<Vector3> stable_vertices, hidden_stable_vertices;
  for (Vertex v: forwardSolver.hullMesh->vertices()){
    if (forwardSolver.vertex_is_stablizable(v)){
      stable_vertices.push_back(normalize(forwardSolver.hullGeometry->inputVertexPositions[v] - forwardSolver.G)+shift);
    }
    else {
      hidden_stable_vertices.push_back(normalize(forwardSolver.hullGeometry->inputVertexPositions[v] - forwardSolver.G)+shift);
    }
  }
  stable_vertices_gm_pc = polyscope::registerPointCloud("stable Vertices Normals", stable_vertices);
  stable_vertices_gm_pc->setPointRadius(pt_cloud_stablizable_radi, false);
  stable_vertices_gm_pc->setPointColor({0.1,0.9,0.1});
  stable_vertices_gm_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

  // hidden stable vertices
  if (show_hidden_stable_vertex_normals){
    hidden_stable_vertices_gm_pc = polyscope::registerPointCloud("hidden stable Vertices Normals", hidden_stable_vertices);
    hidden_stable_vertices_gm_pc->setPointRadius(pt_cloud_stablizable_radi, false);
    hidden_stable_vertices_gm_pc->setPointColor({0.03,0.4,0.03});
    hidden_stable_vertices_gm_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  }
  else {
    if (polyscope::hasPointCloud("hidden stable Vertices Normals"))
      polyscope::removePointCloud("hidden stable Vertices Normals");
  }
}


void draw_stable_face_normals_on_gauss_map(){
  Vector3 shift = {0., gm_distance, 0.};
  std::vector<Vector3> stable_face_normals;
  for (Face f: forwardSolver.hullMesh->faces()){
    Vector3 normal_pos_on_gm = forwardSolver.hullGeometry->faceNormal(f) + shift;
    if (forwardSolver.face_is_stable(f)){
      stable_face_normals.push_back(normal_pos_on_gm);
    }
  }
  stable_face_normals_pc = polyscope::registerPointCloud("stable Face Normals", stable_face_normals);
  stable_face_normals_pc->setPointRadius(face_normal_vertex_gm_radi*1.1, false);
  stable_face_normals_pc->setPointColor({0.9,0.1,0.1});
  stable_face_normals_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere); 
}


void draw_edge_arcs_on_gauss_map(){
  std::vector<std::vector<size_t>> dummy_face{{0,0,0}};
  //    dummy mesh to add curves to
  dummy_psMesh2 = polyscope::registerSurfaceMesh(
      "dummy mesh for gauss map arcs",
      geometry->inputVertexPositions, dummy_face);
  //    add arc per edge
  for (Edge e: forwardSolver.hullMesh->edges()){
    Face f1 = e.halfedge().face(),
         f2 = e.halfedge().twin().face();
    Vector3 n1 = forwardSolver.hullGeometry->faceNormal(f1),
            n2 = forwardSolver.hullGeometry->faceNormal(f2);
    // draw with polyscope
    if (color_arcs) {
      glm::vec3 arc_color = glm::vec3({-1.,0,0}); // default color
      if (forwardSolver.edge_is_stable(e) && forwardSolver.edge_is_stablizable(e)){
        arc_color = both_edge_color;
        draw_arc_on_sphere(n1, n2, gm_shift, gm_radi, arcs_seg_count, e.getIndex(), dummy_psMesh2, 1., arc_color);
      }
      else if (forwardSolver.edge_is_stable(e)){
        arc_color = stable_edge_color;
        draw_arc_on_sphere(n1, n2, gm_shift, gm_radi, arcs_seg_count, e.getIndex(), dummy_psMesh2, 1., arc_color);
      }
      else if (draw_unstable_edge_arcs) {
        if (forwardSolver.edge_is_stablizable(e))
          arc_color = stabilizable_edge_color;
        draw_arc_on_sphere(n1, n2, gm_shift, gm_radi, arcs_seg_count, e.getIndex(), dummy_psMesh2, 1., arc_color);
      }
    } // else; draw all arcs with black
    else {
      draw_arc_on_sphere(n1, n2, gm_shift, gm_radi, arcs_seg_count, e.getIndex(), dummy_psMesh2); // default color
    }
  }
}

void visualize_gauss_map(){
  // just draw the sphere next to the main surface
  std::vector<Vector3> sphere_pos = {gm_shift};
  gauss_map_pc = polyscope::registerPointCloud("Gauss Map", sphere_pos);
  gauss_map_pc->setPointColor({0.74,0.7,0.9});
  gauss_map_pc->setPointRadius(gm_radi, false);
  gauss_map_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

  // point cloud for face normals
  std::vector<Vector3> face_normal_points, stable_face_normals;
  for (Face f: forwardSolver.hullMesh->faces()){
    Vector3 normal_pos_on_gm = forwardSolver.hullGeometry->faceNormal(f) + gm_shift;
    face_normal_points.push_back(normal_pos_on_gm);
    if (forwardSolver.face_is_stable(f)){
      stable_face_normals.push_back(normal_pos_on_gm);
    }
  }
  face_normals_pc = polyscope::registerPointCloud("Face Normals", face_normal_points);
  face_normals_pc->setPointRadius(face_normal_vertex_gm_radi, false);
  face_normals_pc->setPointColor({0.9,0.9,0.9});
  face_normals_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  
  // stable (red) face normals on Gauss map
  draw_stable_face_normals_on_gauss_map();

  // point cloud for stable vertices
  draw_stable_vertices_on_gauss_map();

  // arcs for edge-normals set
  draw_edge_arcs_on_gauss_map();
}


void show_edge_equilibria_on_gauss_map(){
  std::vector<Vector3> edge_equilibria_points, stabilizable_edge_equilibria_points, stable_edge_equilibria_points;
  for (Edge e: forwardSolver.hullMesh->edges()){
    if (forwardSolver.edge_is_stablizable(e) && forwardSolver.edge_is_stable(e)){
      Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
      Vector3 A = forwardSolver.hullGeometry->inputVertexPositions[v1], B = forwardSolver.hullGeometry->inputVertexPositions[v2];
      Vector3 GB = B - G,
              AB = B - A;
      Vector3 ortho_g = GB - AB*dot(AB, GB)/dot(AB,AB);
      edge_equilibria_points.push_back(normalize(ortho_g) + gm_shift);
    }
    else if (forwardSolver.edge_is_stablizable(e)){
      Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
      Vector3 A = forwardSolver.hullGeometry->inputVertexPositions[v1], B = forwardSolver.hullGeometry->inputVertexPositions[v2];
      Vector3 GB = B - G,
              AB = B - A;
      Vector3 ortho_g = GB - AB*dot(AB, GB)/dot(AB,AB);
      stabilizable_edge_equilibria_points.push_back(normalize(ortho_g) + gm_shift);
    }
    else if (forwardSolver.edge_is_stable(e)){
      Vertex v1 = e.firstVertex(), v2 = e.secondVertex();
      Vector3 A = forwardSolver.hullGeometry->inputVertexPositions[v1], B = forwardSolver.hullGeometry->inputVertexPositions[v2];
      Vector3 GB = B - G,
              AB = B - A;
      Vector3 ortho_g = GB - AB*dot(AB, GB)/dot(AB,AB);
      stable_edge_equilibria_points.push_back(normalize(ortho_g) + gm_shift);
    }
  }
  edge_equilibria_pc = polyscope::registerPointCloud("Edge equilibria", edge_equilibria_points);
  edge_equilibria_pc->setPointRadius(face_normal_vertex_gm_radi, false);
  edge_equilibria_pc->setPointColor({0.2, 0.2, 0.9});
  edge_equilibria_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

  if (draw_stable_g_vec_for_unstable_edge_arcs){
    stabilizable_edge_pc = polyscope::registerPointCloud("Stabilizable Edge equilibria", stabilizable_edge_equilibria_points);
    stabilizable_edge_pc->setPointRadius(face_normal_vertex_gm_radi, false);
    stabilizable_edge_pc->setPointColor(stabilizable_edge_color);
    stabilizable_edge_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

    stable_edge_pc = polyscope::registerPointCloud("stable Edge equilibria", stable_edge_equilibria_points);
    stable_edge_pc->setPointRadius(face_normal_vertex_gm_radi, false);
    stable_edge_pc->setPointColor(stable_edge_color);
    stable_edge_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  }
  else{
    if (polyscope::hasPointCloud("stable Edge equilibria"))
      polyscope::removePointCloud("stable Edge equilibria");
    if (polyscope::hasPointCloud("Stabilizable Edge equilibria"))
      polyscope::removePointCloud("Stabilizable Edge equilibria");
  }
}


// visualize center of mass
void draw_G() {
  std::vector<Vector3> G_position = {forwardSolver.G};
  if (polyscope::hasPointCloud("Center of Mass")){
    psG->updatePointPositions(G_position);
  }
  else 
    psG = polyscope::registerPointCloud("Center of Mass", G_position);
  // set some options
  psG->setPointColor({0., 0., 0.});
  psG->setPointRadius(pt_cloud_radi_scale/2.);
  psG->setPointRenderMode(polyscope::PointRenderMode::Sphere);
}

void update_solver(){
  //assuming convex input here
  forwardSolver = Forward3DSolver(mesh, geometry, G);
  // Register the mesh with polyscope
  psInputMesh = polyscope::registerSurfaceMesh(
      "input mesh",
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));
  psInputMesh->setTransparency(0.75);
  draw_G();
}


void visualize_vertex_probabilities(){
  forwardSolver.compute_vertex_probabilities();
  for (Vertex v: forwardSolver.hullMesh->vertices()){
    std::vector<Vector3> positions = {forwardSolver.hullGeometry->inputVertexPositions[v]};
    polyscope::PointCloud* psCloud = polyscope::registerPointCloud("v" + std::to_string(v.getIndex()), positions);
    // set some options
    psCloud->setPointColor({forwardSolver.vertex_probabilities[v], 1., 1.});
    psCloud->setPointRadius(forwardSolver.vertex_probabilities[v] * pt_cloud_radi_scale);
    psCloud->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  }
}


// void visualize_edge_probabilities(){}
void visualize_stable_vertices(){
  std::vector<Vector3> positions;// = {forwardSolver.hullGeometry->inputVertexPositions[v]};
  for (Vertex v: forwardSolver.hullMesh->vertices()){
    if (forwardSolver.vertex_is_stablizable(v)){
      positions.push_back(forwardSolver.hullGeometry->inputVertexPositions[v]);
    }
  }
  polyscope::PointCloud* psCloud = polyscope::registerPointCloud("stabilizable vertices", positions);
  // set some options
  psCloud->setPointColor({0.1, .9, .1});
  psCloud->setPointRadius(pt_cloud_stablizable_radi);
  psCloud->setPointRenderMode(polyscope::PointRenderMode::Sphere);
}

void visualize_edge_stability(){
  std::vector<std::vector<size_t>> dummy_face{{1,1,1}};
  dummy_psMesh1 = polyscope::registerSurfaceMesh(
      "dummy mesh for edges",
      geometry->inputVertexPositions, dummy_face);

  std::vector<std::array<size_t, 2>> stable_edgeInds, stablilizable_edgeInds, both_edgeInds;
  std::vector<Vector3> stable_positions, stablilizable_positions, both_positions;
  size_t stable_counter = 0, stablizable_counter = 0, both_counter = 0;
  for (Edge e: forwardSolver.hullMesh->edges()){
    Vector3 p1 = forwardSolver.hullGeometry->inputVertexPositions[e.firstVertex()],
            p2 = forwardSolver.hullGeometry->inputVertexPositions[e.secondVertex()];
    size_t flag = 0;
    if (forwardSolver.edge_is_stable(e)){
      stable_positions.push_back(p1); stable_positions.push_back(p2);
      stable_edgeInds.push_back({stable_counter, stable_counter + 1});
      stable_counter += 2;
      flag++;
    }
    if (forwardSolver.edge_is_stablizable(e)){
      stablilizable_positions.push_back(p1); stablilizable_positions.push_back(p2);
      stablilizable_edgeInds.push_back({stablizable_counter, stablizable_counter + 1});
      stablizable_counter += 2;
      flag++;
    }
    if (flag == 2){
      both_positions.push_back(p1); both_positions.push_back(p2);
      both_edgeInds.push_back({both_counter, both_counter + 1});
      both_counter += 2;
    }
  }
  polyscope::SurfaceGraphQuantity* psStableEdges =  dummy_psMesh1->addSurfaceGraphQuantity("stable Edges", stable_positions, stable_edgeInds);
  psStableEdges->setRadius(stable_edge_radi, true);
  psStableEdges->setColor(stable_edge_color);
  psStableEdges->setEnabled(true);
  polyscope::SurfaceGraphQuantity* psStablilizableEdges =  dummy_psMesh1->addSurfaceGraphQuantity("stablizable Edges", stablilizable_positions, stablilizable_edgeInds);
  psStablilizableEdges->setRadius(stablizable_edge_radi, true);
  psStablilizableEdges->setColor(stabilizable_edge_color);
  psStablilizableEdges->setEnabled(true);
  polyscope::SurfaceGraphQuantity* psBothEdges =  dummy_psMesh1->addSurfaceGraphQuantity("stable && stablizable Edges", both_positions, both_edgeInds);
  psBothEdges->setRadius(both_edge_radi, true);
  psBothEdges->setColor(both_edge_color);
  psBothEdges->setEnabled(true);
}


void visualize_face_stability(){
  std::vector<std::array<double, 3>> fColor(forwardSolver.hullMesh->nFaces());
  for (Face f: forwardSolver.hullMesh->faces()){
    if (forwardSolver.face_is_stable(f))
      fColor[f.getIndex()] = {1., 0.1, 0.1};
    else
      fColor[f.getIndex()] = {0.9, 0.9, 1.};
  }
  polyscope::SurfaceFaceColorQuantity *faceQnty = psInputMesh->addFaceColorQuantity("face stability", fColor);
  faceQnty->setEnabled(true);
}


// generate simple examples
void generate_polyhedron_example(std::string poly_str){
    // readManifoldSurfaceMesh()
    std::tie(mesh_ptr, geometry_ptr) = generate_polyhedra(poly_str);
    mesh = mesh_ptr.release();
    geometry = geometry_ptr.release();
}

// initialize Markov chain and do the pre-computes
void initiate_markov_model(){
  markov_model = new RollingMarkovModel(&forwardSolver);
  markov_model->initialize_pre_computes();
  markov_model->split_chain_edges();
}


// visualize 
void visualize_sudo_faces(){
  std::vector<Vector3> sf_positions;
  size_t sf_count = 0;
  for (Halfedge he: markov_model->mesh->halfedges()){
    SudoFace *curr_sf = markov_model->root_sudo_face[he];
    if (curr_sf != nullptr){
      while (curr_sf->next_sudo_face != curr_sf){
        sf_positions.push_back(curr_sf->normal + gm_shift);
        sf_count++;
        curr_sf = curr_sf->next_sudo_face;
      }
    }
  }
  
  sudo_faces_pc = polyscope::registerPointCloud("sudo faces", sf_positions);
  sudo_faces_pc->setEnabled(true);
  sudo_faces_pc->setPointRadius(face_normal_vertex_gm_radi * 1.2, false);
  sudo_faces_pc->setPointColor({0.,1.,1.});
}


// color input polyhedra with default
void color_faces_with_default(){
  FaceData<Vector3> face_colors(*forwardSolver.hullMesh, default_face_color);
  polyscope::SurfaceFaceColorQuantity *fColor = psInputMesh->addFaceColorQuantity("state vis color", face_colors);
  fColor->setEnabled(true);
}


// show current g vector
void visualize_g_vec(){
  std::vector<Vector3> the_g_vec = {forwardSolver.curr_g_vec};
  polyscope::PointCloudVectorQuantity *psG_vec = psG->addVectorQuantity("g_vec", the_g_vec);
  psG_vec->setEnabled(true);
  psG_vec->setVectorRadius(curve_radi_scale * 1.);
  psG_vec->setVectorLengthScale(0.2);
  psG_vec->setVectorColor({0.8,0.1,0.1});
}


// visualize current touching element: Vertex/Edge/Face
void visualize_contact(){
  if (polyscope::hasPointCloud("current Vertex")) polyscope::removePointCloud("current Vertex");
  dummy_forward_vis->removeQuantity("current contact edge"); // has a built-in existance checker
  if (forwardSolver.curr_f.getIndex() == INVALID_IND) color_faces_with_default();
  // add the other two
  if (forwardSolver.curr_v.getIndex() != INVALID_IND) {
    printf("at Vertex\n");
    // show contact vertex on polyhedra
    std::vector<Vector3> curr_state_pos = {forwardSolver.hullGeometry->inputVertexPositions[forwardSolver.curr_v]}; // first and second should be the same since we just initialized.
    curr_state_pt = polyscope::registerPointCloud("current Vertex", curr_state_pos);
    curr_state_pt->setEnabled(true);
    curr_state_pt->setPointRadius(pt_cloud_radi_scale/2.);
    
    // show gravity vec on Gauss Map
    std::vector<Vector3> curr_g_vec_pos = {forwardSolver.curr_g_vec + gm_shift}; // first and second should be the same since we just initialized.
    curr_g_vec_gm_pt = polyscope::registerPointCloud("current g vec", curr_g_vec_pos);
    curr_g_vec_gm_pt->setEnabled(true);
    curr_g_vec_gm_pt->setPointRadius(face_normal_vertex_gm_radi * 1.1, false);
    curr_g_vec_gm_pt->setPointColor({0.,0.,0.});
  }
  else if (forwardSolver.curr_e.getIndex() != INVALID_IND){
    printf("at Edge\n");
    Vertex v1 = forwardSolver.curr_e.firstVertex(),
           v2 = forwardSolver.curr_e.secondVertex();
    Vector3 p1 = forwardSolver.hullGeometry->inputVertexPositions[v1],
            p2 = forwardSolver.hullGeometry->inputVertexPositions[v2];
    std::vector<std::array<size_t, 2>> edgeInds;
    std::vector<Vector3> positions;
    edgeInds.push_back({0, 1});
    positions.push_back(p1); positions.push_back(p2);
    curr_state_segment =  dummy_forward_vis->addSurfaceGraphQuantity("current contact edge", positions, edgeInds);
    curr_state_segment->setRadius(curve_radi_scale/2.3);
    curr_state_segment->setColor({0., 0., 1.});
    curr_state_segment->setEnabled(true);
  }
  else if (forwardSolver.curr_f.getIndex() != INVALID_IND){
    face_colors = FaceData<Vector3>(*forwardSolver.hullMesh, default_face_color);
    face_colors[forwardSolver.curr_f] = curr_face_color;
    polyscope::SurfaceFaceColorQuantity *fColor = psInputMesh->addFaceColorQuantity("state vis color", face_colors);
    fColor->setEnabled(true);
  }
  else {
    polyscope::warning("all elements are Invalid??\n");
  }
  // TODO: maybe add Snail trail?
}


void initialize_state_vis(){
  // for later single segment curve addition (for stable edges)
  std::vector<std::vector<size_t>> dummy_face{{0,0,0}};
  std::vector<Vector3> dummy_pos = {Vector3({0.,0.,0.})};
  dummy_forward_vis = polyscope::registerSurfaceMesh("state check mesh", dummy_pos, dummy_face); // nothing matters in this line
  // will add curves to this later
  draw_G();
  visualize_contact();
  visualize_g_vec();
  color_faces_with_default();
}


void color_faces(){
  if (recolor_faces){
    face_colors = generate_random_colors(forwardSolver.hullMesh);
    for (Face f: forwardSolver.hullMesh->faces()){
      if (!forwardSolver.face_is_stable(f))
        face_colors[f] = default_face_color;
    }
    recolor_faces = false;
  }
}


void visualize_colored_polyhedra(){
  VertexData<Vector3> shifted_positions(*forwardSolver.hullMesh);
  for (Vertex v: forwardSolver.hullMesh->vertices()){
    shifted_positions[v] = forwardSolver.hullGeometry->inputVertexPositions[v] + colored_shift;
  }
  coloredPsMesh = polyscope::registerSurfaceMesh("colored polyhedra", shifted_positions, forwardSolver.hullMesh->getFaceVertexList());
  // generate random colors and color the faces
  polyscope::SurfaceFaceColorQuantity *faceQnty = coloredPsMesh->addFaceColorQuantity("random face colors", face_colors);
  faceQnty->setEnabled(true);
  // // add colors to the original polyhedra as well?
  // polyscope::SurfaceFaceColorQuantity *faceQnty2 = psMesh->addFaceColorQuantity("random face colors2", face_colors);
  // faceQnty2->setEnabled(true);
}


void update_visuals_with_G(){
  forwardSolver.G = G;
  draw_G();
  // stuff on the polyhedra
  visualize_edge_stability();
  visualize_face_stability();
  visualize_stable_vertices();
  
  // coloring stuff
  recolor_faces = true;// so bad lol
  color_faces();
  visualize_colored_polyhedra();
  // Gauss map stuff
  if(color_arcs){
    draw_edge_arcs_on_gauss_map();
  }
  draw_stable_vertices_on_gauss_map();
  show_edge_equilibria_on_gauss_map();
  draw_stable_face_normals_on_gauss_map();
  if (polyscope::hasSurfaceMesh("dummy mesh for gauss map arcs") && draw_stable_patches)
    draw_stable_patches_on_gauss_map();
  initiate_markov_model();
  visualize_sudo_faces();
}


// maybe move to another file, to merge with bullet sim; this and some other functions
// sample and raster; colors should be generated beforehand
void build_raster_image(){
  FaceData<std::vector<Vector3>> face_samples(*forwardSolver.hullMesh);
  int total_invalids = 0;
  // need to re-iterate later with the same order
  FaceData<std::vector<Vector3>> initial_rolling_dirs(*forwardSolver.hullMesh);
  for (int i = 0; i < sample_count; i++){
    Vector3 random_g_vec = {randomReal(-1,1), randomReal(-1,1), randomReal(-1,1)};
    if (random_g_vec.norm() <= 1){
      if (i % 5000 == 0)
        printf("$$$ at sample %d\n", i);
      random_g_vec /= norm(random_g_vec);
      Face touching_face = forwardSolver.final_touching_face(random_g_vec);
      if (touching_face.getIndex() == INVALID_IND){
        total_invalids++;
        continue;
      }
      face_samples[touching_face].push_back(random_g_vec);
      if (normalize_vecF)
        initial_rolling_dirs[touching_face].push_back(forwardSolver.initial_roll_dir.normalize());
      else 
        initial_rolling_dirs[touching_face].push_back(forwardSolver.initial_roll_dir);
    }
  }
  printf(" ###### total invalid faces: %d  ######\n", total_invalids);
  std::vector<Vector3> raster_positions,
                       raster_colors;
  std::vector<Vector3> all_rolling_dirs;
  for (Face f: forwardSolver.hullMesh->faces()){
    std::vector<Vector3> tmp_points = face_samples[f];
    for (Vector3 tmp_p: tmp_points){
      raster_positions.push_back(tmp_p + gm_shift);
      raster_colors.push_back(face_colors[f]);
    }
    std::vector<Vector3> tmp_rolling_dirs = initial_rolling_dirs[f];
    for (Vector3 rolling_dir: tmp_rolling_dirs){
      all_rolling_dirs.push_back(rolling_dir);
    }
  }
  raster_pc = polyscope::registerPointCloud("raster point cloud", raster_positions);
  raster_pc->setPointRadius(0.0078, false);
  polyscope::PointCloudColorQuantity* pc_col_quant = raster_pc->addColorQuantity("random color", raster_colors);
  pc_col_quant->setEnabled(true);
  polyscope::PointCloudVectorQuantity* raster_pc_vec_field = raster_pc->addVectorQuantity("init roll directions", all_rolling_dirs);
  raster_pc_vec_field->setEnabled(true);
  raster_pc_vec_field->setVectorLengthScale(0.05, false);
  raster_pc_vec_field->setVectorRadius(0.003, false);
  raster_pc_vec_field->setVectorColor(vecF_color);
}



// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
  if (ImGui::BeginCombo("##combo1", all_polygons_current_item.c_str())){ // The second parameter is the label previewed before opening the combo.
      for (std::string tmp_str: all_polyhedra_items){ // This enables not having to have a const char* arr[]. Or maybe I'm just a noob.
          bool is_selected = (all_polygons_current_item == tmp_str.c_str()); // You can store your selection however you want, outside or inside your objects
          if (ImGui::Selectable(tmp_str.c_str(), is_selected)){ // selected smth
              all_polygons_current_item = tmp_str;
              generate_polyhedron_example(all_polygons_current_item);
              update_solver();
              recolor_faces = true;
              color_faces();
          }
          if (is_selected)
              ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
      }
      ImGui::EndCombo();
  }
  if (ImGui::Button("visualize vertex probabilities"))
    visualize_vertex_probabilities();
  
  if (ImGui::SliderFloat("vertex radi scale", &pt_cloud_radi_scale, 0., 1.)){
    draw_G();
  }

  if (ImGui::Button("show edge stability status")  ||
      ImGui::SliderFloat("stable curve radi", &stable_edge_radi, 0., 0.05)||
      ImGui::SliderFloat("stablizable curve radi", &stablizable_edge_radi, 0., 0.05) ||
      ImGui::SliderFloat("both curve radi", &both_edge_radi, 0., 0.05)){
      visualize_edge_stability();
      visualize_face_stability();
  }
  if (ImGui::Button("show face stability status")){
      visualize_face_stability();
  }
  if (ImGui::Button("show vertex stability status") || 
      ImGui::SliderFloat("stable vertices radi", &pt_cloud_stablizable_radi, 0., 0.05)){
      visualize_stable_vertices();
  }
  
  // if (ImGui::SliderFloat("edge radi scale", &curve_radi_scale, 0., 1.)) 
  //   visualize_edge_probabilities();
  
  // if (ImGui::Button("visualize edge probabilities"))
  //   visualize_edge_probabilities();

  if (ImGui::SliderFloat("G radi scale", &G_r, -3., 3.)||
      ImGui::SliderFloat("G theta", &G_theta, 0., 2*PI)||
      ImGui::SliderFloat("G phi", &G_phi, 0., 2*PI)) {
    G = {cos(G_phi)*sin(G_theta)*G_r, cos(G_phi)*cos(G_theta)*G_r, sin(G_phi)*G_r};
    if (G_is_inside(*forwardSolver.hullMesh, *forwardSolver.hullGeometry, G)){
      update_visuals_with_G();
      if (real_time_raster){
        color_faces();
        build_raster_image();
      }
    }
  }
  if (ImGui::Button("uniform mass G")){
    G = find_center_of_mass(*mesh, *geometry);
    update_visuals_with_G();
  }
  if (ImGui::Button("Show Gauss Map") 
      // || 
      // ImGui::SliderInt("seg count for arcs", &arcs_seg_count, 1, 100)||
      // ImGui::SliderFloat("arc curve radi", &arc_curve_radi, 0., 0.04)||
      // ImGui::SliderFloat("face normal vertex radi", &face_normal_vertex_gm_radi, 0., 0.04)
      ){///face_normal_vertex_gm_radi
    visualize_gauss_map();
  }
  if (ImGui::Checkbox("colored arcs (slows)", &color_arcs));
  if (ImGui::Checkbox("draw unstable edge arcs", &draw_unstable_edge_arcs)) draw_edge_arcs_on_gauss_map();
  if (ImGui::Checkbox("show to-be stable g_vec for unstable edge arcs", &draw_stable_g_vec_for_unstable_edge_arcs)) show_edge_equilibria_on_gauss_map();
  if (ImGui::Checkbox("show hidden stable vertex normals", &show_hidden_stable_vertex_normals)) draw_stable_vertices_on_gauss_map();
  if (ImGui::Checkbox("draw stable patches", &draw_stable_patches)) draw_stable_patches_on_gauss_map();
  if (ImGui::Button("Draw patches")){
    draw_stable_patches_on_gauss_map();
  }
  
  // my simulation 
  if (ImGui::Button("initialize g_vec") ||
    ImGui::SliderFloat("initial g_vec theta", &g_vec_theta, 0., 2*PI)||
    ImGui::SliderFloat("initial g_vec phi", &g_vec_phi, 0., 2*PI)) {
    initial_g_vec = {cos(g_vec_phi)*sin(g_vec_theta), cos(g_vec_phi)*cos(g_vec_theta), sin(g_vec_phi)};
    forwardSolver.find_contact(initial_g_vec);
    initialize_state_vis();
    // snail trail stuff
    snail_trail_dummy_counter = 0;
    std::vector<std::vector<size_t>> dummy_face{{0,0,0}};
    dummy_ps_mesh_for_snail_trail = polyscope::registerSurfaceMesh("dummy mesh snail trail", geometry->inputVertexPositions, dummy_face);
  }
  if (ImGui::Button("next state")){
    old_g_vec = forwardSolver.curr_g_vec;
    forwardSolver.next_state();
    new_g_vec = forwardSolver.curr_g_vec;
    visualize_g_vec();
    visualize_contact();
    if(show_snail_trail && norm(old_g_vec-new_g_vec) != 0.){ // proly dont have to use tol
      draw_arc_on_sphere(old_g_vec, new_g_vec, gm_shift, gm_radi, arcs_seg_count, 200 + snail_trail_dummy_counter, dummy_ps_mesh_for_snail_trail,
                         1.5);
      snail_trail_dummy_counter++;
    }
  }
  if (ImGui::SliderInt("sample count", &sample_count, 1000, 100000));
  if (ImGui::Button("build raster image")){
    color_faces();
    visualize_colored_polyhedra();
    build_raster_image();
  }
  if(ImGui::Button("recolor faces")){
    recolor_faces = true;// so bad lol
    color_faces();
    visualize_colored_polyhedra();
  }
  if (ImGui::Checkbox("real time raster", &real_time_raster));
  if (ImGui::Checkbox("normalize vector field", &normalize_vecF));

  if (ImGui::Button("show sudo faces")){
    initiate_markov_model();
    visualize_sudo_faces();
  }
}


int main(int argc, char **argv) {
  // build mesh
  generate_polyhedron_example(all_polygons_current_item);
  G = {0.,0,0};
  update_solver();
  // build the solver
  

  // Initialize polyscope
  polyscope::init();
  color_faces();
  visualize_colored_polyhedra();
  // Set the callback function
  polyscope::state::userCallback = myCallback;
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  polyscope::view::upDir = polyscope::view::UpDir::NegZUp;
  
  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
