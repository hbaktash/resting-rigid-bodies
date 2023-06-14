#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
// #include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include "args/args.hxx"
#include "imgui.h"

#include "forward3D.h"

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
        initial_g_vec;

Forward3DSolver forwardSolver;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psInputMesh, *dummy_psMesh1, 
                       *dummy_psMesh2, *dummy_psMesh3, *dummy_forward_vis;

polyscope::PointCloud *psG, *curr_state_pt, 
                      *gauss_map_pc, *face_normals_pc; // point cloud with single G

polyscope::SurfaceGraphQuantity* curr_state_segment;

float pt_cloud_radi_scale = 0.1,
      curve_radi_scale = 0.1,
      G_r = 0.,
      G_theta = 0.,
      G_phi = 0.;

float stable_edge_radi = 0.0,
      stablizable_edge_radi = 0.03,
      both_edge_radi = 0.02,
      pt_cloud_stablizable_radi = 0.03,
      arc_curve_radi = 0.01,
      face_normal_vertex_gm_radi = 0.03;

double gm_distance = 2.1,
       gm_radi = 1.0;

int sample_count = 1000,
    arcs_seg_count = 13;

// example choice
std::vector<std::string> all_polyhedra_items = {std::string("cube"), std::string("tet"), std::string("sliced tet")};
std::string all_polygons_current_item = "tet";
static const char* all_polygons_current_item_c_str = "tet";


Vector3 spherical_to_xyz(double r, double phi, double theta){
  return Vector3({r*cos(phi)*sin(theta), r*cos(phi)*cos(theta), r*sin(phi)});
}

// draw an arc connecting two points on the sphere; for Gauss map purposes
void draw_arc_on_sphere(Vector3 p1, Vector3 p2, Vector3 center, double radius, size_t seg_count, size_t edge_ind){
  // p1, p2 just represent normal vectors
  if (norm(p1) > 1.01 || norm(p2) > 1.01)
    polyscope::warning("wtf? p1, p2 norm larger than 1");
  
  std::vector<std::array<size_t, 2>> edgeInds;
  std::vector<Vector3> positions;
  double sqrt_radi = sqrt(radius);
  // walk on p1-p2 segment
  Vector3 curr_point = p1,
          forward_vec = (p2-p1)/(double)seg_count;
  Vector3 next_point = curr_point + forward_vec;
  Vector3 curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center ,
          next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
  positions.push_back(curr_point_on_sphere);
  for (size_t i = 0; i < seg_count; i++){
    // add to positions list
    curr_point_on_sphere = normalize(curr_point) * sqrt_radi + center ,
    next_point_on_sphere = normalize(next_point) * sqrt_radi + center;
    positions.push_back(next_point_on_sphere);
    // add segment indices
    edgeInds.push_back({i, i+1});

    // update points
    curr_point = next_point;
    next_point += forward_vec;
  }
  polyscope::SurfaceGraphQuantity* psArcCurve = dummy_psMesh2->addSurfaceGraphQuantity("Arc curve " + std::to_string(edge_ind), positions, edgeInds);
  psArcCurve->setRadius(arc_curve_radi, false);
  psArcCurve->setColor({0.03, 0.03, 0.03});
  psArcCurve->setEnabled(true);
}


void visualize_gauss_map(){
  // just draw the sphere next to the main surface
  Vector3 shift = {0., gm_distance, 0.};
  std::vector<Vector3> sphere_pos = {shift};
  gauss_map_pc = polyscope::registerPointCloud("Gauss Map", sphere_pos);
  gauss_map_pc->setPointColor({0.74,0.7,0.9});
  gauss_map_pc->setPointRadius(gm_radi, false);
  gauss_map_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);

  // point cloud for face normals
  std::vector<Vector3> face_normal_points;
  for (Face f: forwardSolver.hullMesh->faces()){
    face_normal_points.push_back(forwardSolver.hullGeometry->faceNormal(f) + shift);
  }
  face_normals_pc = polyscope::registerPointCloud("Face Normals", face_normal_points);
  face_normals_pc->setPointRadius(face_normal_vertex_gm_radi, false);
  face_normals_pc->setPointColor({0.,0.2,0.99});
  face_normals_pc->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  
  // arcs for edge normals set
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
    draw_arc_on_sphere(n1, n2, shift, gm_radi, arcs_seg_count, e.getIndex());
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
  psG->setPointColor({1., 1., 1.});
  psG->setPointRadius(pt_cloud_radi_scale/2.);
  psG->setPointRenderMode(polyscope::PointRenderMode::Sphere);
}

void update_solver(){
  forwardSolver = Forward3DSolver(mesh, geometry, G);
  //assuming convex input here
  forwardSolver.hullMesh = new ManifoldSurfaceMesh(mesh->getFaceVertexList()); //mesh->copy().release();
  forwardSolver.hullGeometry = new VertexPositionGeometry(*forwardSolver.hullMesh); // geometry->copy().release();
  for (Vertex v: forwardSolver.hullMesh->vertices()){
    forwardSolver.hullGeometry->inputVertexPositions[v] = geometry->inputVertexPositions[v.getIndex()];
  }
  // forwardSolver.compute_vertex_probabilities();

  // Register the mesh with polyscope
  psInputMesh = polyscope::registerSurfaceMesh(
      "input mesh",
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));
  psInputMesh->setTransparency(0.95);
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


void visualize_vec_from_G(Vector3 vec){
  std::vector<Vector3> the_g_vec = {vec};
  polyscope::PointCloudVectorQuantity *psG_vec = psG->addVectorQuantity("g_vec", the_g_vec);
  psG_vec->setEnabled(true);
  psG_vec->setVectorRadius(curve_radi_scale * 1.);
  psG_vec->setVectorLengthScale((dot(vec, vec)));
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
  psStableEdges->setColor({0., 1., 1.});
  psStableEdges->setEnabled(true);
  polyscope::SurfaceGraphQuantity* psStablilizableEdges =  dummy_psMesh1->addSurfaceGraphQuantity("stablizable Edges", stablilizable_positions, stablilizable_edgeInds);
  psStablilizableEdges->setRadius(stablizable_edge_radi, true);
  psStablilizableEdges->setColor({0.1, 0.9, 0.2});
  psStablilizableEdges->setEnabled(true);
  polyscope::SurfaceGraphQuantity* psBothEdges =  dummy_psMesh1->addSurfaceGraphQuantity("stable && stablizable Edges", both_positions, both_edgeInds);
  psBothEdges->setRadius(both_edge_radi, true);
  psBothEdges->setColor({0.2, 0.2, 0.9});
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
    std::vector<std::vector<size_t>> faces;
    std::vector<Vector3> positions;
    int n;
    if (std::strcmp(poly_str.c_str(), "tet") == 0){
      n = 4;
      faces = {{0, 2, 1},
                              {0, 1, 3},
                              {0, 3, 2},
                              {2, 3, 1}};
      double theta0 = 0., theta1 = 2.*PI/3., theta2 = 4.*PI/3.,
             phi0 = PI/6.,
             theta3 = 0., 
             phi1 = -PI/2.;
      
      positions = { spherical_to_xyz(1., phi0, theta0),
                    spherical_to_xyz(1., phi0, theta1),
                    spherical_to_xyz(1., phi0, theta2),
                    spherical_to_xyz(1., phi1, theta3)};
    }
    else if (std::strcmp(poly_str.c_str(), "sliced tet") == 0){
      double theta0 = 0., theta1 = 2.*PI/3., theta2 = 4.*PI/3.,
              phi0 = PI/6.,
              theta3 = 0., 
              phi1 = -PI/2.;
      Vector3 p1 = spherical_to_xyz(1., phi0, theta0), //0
              p2 = spherical_to_xyz(1., phi0, theta1), // gone
              p3 = spherical_to_xyz(1., phi0, theta2), // gone
              p = spherical_to_xyz(1., phi1, theta3); // 1
      Vector3 p2p  = 0.8*p2 + 0.2*p, // 2
              p2p3 = 0.7*p2 + 0.3*p3, // 3
              p2p1 = 0.8*p2 + 0.2*p1; // 4
      Vector3 p3p  = 0.8*p3 + 0.2*p, // 5
              p3p2 = 0.7*p3 + 0.3*p2, // 6
              p3p1 = 0.8*p3 + 0.2*p1; // 7
      // making the quad faces non-planar
      double shrinkage = 0.7;
      p2p *= 0.8; 
      p2p1 *= shrinkage;
      p3p *= 0.8; 
      p3p1 *= shrinkage;
      positions.push_back(p); positions.push_back(p1);
      positions.push_back(p2p); positions.push_back(p2p3); positions.push_back(p2p1); 
      positions.push_back(p3p); positions.push_back(p3p2); positions.push_back(p3p1);
      faces = {{0, 3, 6},
               {0, 2, 3},
               {0, 6, 5},
               {0, 1, 2},
               {0, 5, 1},
               {1, 3, 4},
               {1, 7, 6},   
               {1, 6, 3},
               {5, 6, 7},
               {2, 4, 3},
               {1, 4, 2},
               {1, 5, 7}};
    }
    else if (std::strcmp(poly_str.c_str(), "rndbs 6-gon 1") == 0){
    }
    else if (std::strcmp(poly_str.c_str(), "rndbs 9-gon 1") == 0){
    }
    else {
      throw std::runtime_error("no valid string provided\n");
    }
    // std::unique_ptr<ManifoldSurfaceMesh> poly_triangulated;
    // poly_triangulated.reset(new ManifoldSurfaceMesh(greedy_triangulation));
    mesh = new ManifoldSurfaceMesh(faces);
    geometry = new VertexPositionGeometry(*mesh);
    for (Vertex v : mesh->vertices()) {
        // Use the low-level indexers here since we're constructing
        // printf("v %d\n", v.getIndex());
        geometry->inputVertexPositions[v] = positions[v.getIndex()];
    }
}


void visualize_contact_point(){
  // single point cloud for the single contact point
  if (forwardSolver.curr_state.first == forwardSolver.curr_state.second){
    std::vector<Vector3> curr_state_pos = {forwardSolver.hullGeometry->inputVertexPositions[forwardSolver.curr_state.first]}; // first and second should be the same since we just initialized.
    curr_state_pt = polyscope::registerPointCloud("current state", curr_state_pos);
    curr_state_pt->setEnabled(true);
    curr_state_pt->setPointRadius(pt_cloud_radi_scale/2.);
  }
}

void initialize_state_vis(){
  visualize_vec_from_G(forwardSolver.current_g_vec);
  visualize_contact_point();
  draw_G();
  // for later single segment curve addition
  std::vector<std::vector<size_t>> dummy_face{{0,0,0}};
  std::vector<Vector3> dummy_pos = {Vector3({0.,0.,0.})};
  dummy_forward_vis = polyscope::registerSurfaceMesh("state check mesh", dummy_pos, dummy_face); // nothing matters in this line
  // will add curves to this later
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
    forwardSolver.G = G;
    draw_G();
    visualize_edge_stability();
    visualize_face_stability();
    visualize_stable_vertices();
  }

  if (ImGui::Button("Show Gauss Map") || 
      ImGui::SliderInt("seg count for arcs", &arcs_seg_count, 1, 100)||
      ImGui::SliderFloat("arc curve radi", &arc_curve_radi, 0., 0.04)||
      ImGui::SliderFloat("face normal vertex radi", &face_normal_vertex_gm_radi, 0., 0.04)){///face_normal_vertex_gm_radi
    visualize_gauss_map();
  }
    
  // if (ImGui::InputInt("compute empirical edge probabilities", &sample_count, 100)) {
  //   forwardSolver.build_next_edge_tracer();
  //   forwardSolver.empirically_build_probabilities(sample_count);
  // }
  // if (ImGui::SliderFloat("g-vec angle", &g_vec_angle, 0., 2*PI)) {
    // initial_g_vec = {cos(g_vec_angle), sin(g_vec_angle), 0};
    // forwardSolver.find_contact(initial_g_vec);
    // forwardSolver.build_next_edge_tracer();
    // initialize_state_vis();
  // }
  // if (ImGui::Button("reset state")) {
    // forwardSolver.find_contact(initial_g_vec);
    // forwardSolver.build_next_edge_tracer();
    // initialize_state_vis();
  // }
  // if (ImGui::Button("next state")) {
    // forwardSolver.next_state();
    // visualize_g_vec();
    // visualize_contact_point();
  // }
}


int main(int argc, char **argv) {

  // Configure the argument parser
  // args::ArgumentParser parser("geometry-central & Polyscope example project");
  // args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // // Parse args
  // try {
  //   parser.ParseCLI(argc, argv);
  // } catch (args::Help &h) {
  //   std::cout << parser;
  //   return 0;
  // } catch (args::ParseError &e) {
  //   std::cerr << e.what() << std::endl;
  //   std::cerr << parser;
  //   return 1;
  // }

  // // Make sure a mesh name was given
  // if (!inputFilename) {
  //   std::cerr << "Please specify a mesh file as argument" << std::endl;
  //   return EXIT_FAILURE;
  // }

  // build mesh
  generate_polyhedron_example(all_polygons_current_item);
  G = {0.,0,0};
  update_solver();
  // build the solver
  
  // Load mesh
  // std::tie(mesh, geometry) = readManifoldSurfaceMesh(args::get(inputFilename));

  // Initialize polyscope
  polyscope::init();
  
  // Set the callback function
  polyscope::state::userCallback = myCallback;
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  polyscope::view::upDir = polyscope::view::UpDir::NegZUp;

  // Set vertex tangent spaces
  
  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
