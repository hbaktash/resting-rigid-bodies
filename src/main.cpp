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
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include "args/args.hxx"
#include "imgui.h"

#include "forward.h"

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

ForwardSolver forwardSolver;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psInputMesh, *dummy_psMesh1, *dummy_psMesh2, *dummy_psMesh3, 
                       *dummy_forward_vis;

polyscope::PointCloud *psG, *curr_state_pt; // point cloud with single G


float pt_cloud_radi_scale = 0.1,
      curve_radi_scale = 0.1,
      G_scale = 1.,
      G_angle = 0.,
      g_vec_angle = 0.;
int sample_count = 1e3;

// example choice
std::vector<std::string> all_polygon_items = {std::string("cube"), std::string("skewed 4-gon"), std::string("rndbs 6-gon 1"), std::string("rndbs 9-gon 1")};
std::string all_polygons_current_item = "cube";
static const char* all_polygons_current_item_c_str = "cube";


void visualize_boundary_curves(){
  std::vector<std::vector<size_t>> dummy_face{{1,1,1}};
  auto positions = forwardSolver.hullGeometry->inputVertexPositions.toVector();
  std::vector<std::array<size_t, 2>> edgeInds;
  for (Edge e: forwardSolver.hullMesh->edges()){
    if (e.isBoundary()){
      edgeInds.push_back({e.firstVertex().getIndex(), e.secondVertex().getIndex()});
    }
  }
  auto convex_bnd_net = polyscope::registerCurveNetwork("convex boundary", positions, edgeInds);
}


void visualize_vertex_probabilities(){
  
  for (Vertex v: forwardSolver.hullMesh->vertices()){
    std::vector<Vector3> positions = {forwardSolver.hullGeometry->inputVertexPositions[v]};
    polyscope::PointCloud* psCloud = polyscope::registerPointCloud("v" + std::to_string(v.getIndex()), positions);
    // set some options
    psCloud->setPointColor({forwardSolver.vertex_probabilities[v], 1., 1.});
    psCloud->setPointRadius(forwardSolver.vertex_probabilities[v] * pt_cloud_radi_scale);
    psCloud->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  }

}

void visualize_edge_probabilities(){
  forwardSolver.build_next_edge_tracer();
  forwardSolver.compute_initial_edge_probabilities();
  forwardSolver.compute_final_edge_probabilities();

  // initial probs
  // auto positions = forwardSolver.hullGeometry->inputVertexPositions.toVector();
  // std::vector<std::array<size_t, 2>> edgeInds;
  for (Edge e: forwardSolver.hullMesh->edges()){
    if (e.isBoundary()){
      std::vector<std::array<size_t, 2>> edgeInds;
      std::vector<Vector3> positions;
      edgeInds.push_back({0, 1});
      Vector3 p1 = forwardSolver.hullGeometry->inputVertexPositions[e.firstVertex()],
              p2 = forwardSolver.hullGeometry->inputVertexPositions[e.secondVertex()];
      positions.push_back(p1); positions.push_back(p2);
      // initial prob visualize
      auto psSegmentNet = polyscope::registerCurveNetwork("initial prob e" + std::to_string(e.getIndex()), positions, edgeInds);
      psSegmentNet->setRadius(forwardSolver.initial_edge_probabilities[e]*curve_radi_scale);
      psSegmentNet->setColor({1., forwardSolver.initial_edge_probabilities[e], 1.});
      psSegmentNet->setEnabled(true);
      // final prob visualize
      auto psSegment2Net =  polyscope::registerCurveNetwork("final prob e" + std::to_string(e.getIndex()), positions, edgeInds);
      psSegment2Net->setRadius(forwardSolver.final_edge_probabilities[e]*curve_radi_scale);
      psSegment2Net->setColor({1., 1., forwardSolver.final_edge_probabilities[e]});
      psSegment2Net->setEnabled(true);
    }
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

// generate simple examples
void generate_polygon_example(std::string poly_str){
    std::vector<std::vector<size_t>> greedy_triangulation;
    std::vector<Vector3> positions;
    int n;
    if (std::strcmp(poly_str.c_str(), "cube") == 0){
      n = 4;
      greedy_triangulation = {{0, 1, 2},
                              {0, 2, 3}};
      positions = {{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}};
    }
    else if (std::strcmp(poly_str.c_str(), "skewed 4-gon") == 0){
      n = 4;
      greedy_triangulation = {{0, 1, 2},
                              {0, 2, 3}};
      positions = {{-1*3., -1*3., 0}, {1*1.5, -1*1.5, 0}, {1*2, 1*2, 0}, {-1, 1, 0}};
    }
    else if (std::strcmp(poly_str.c_str(), "rndbs 6-gon 1") == 0){
      n = 6;
      greedy_triangulation = {{0, 1, 2},
                              {0, 2, 3},
                              {0, 3, 4},
                              {0, 4, 5}};
      for (int i = 0; i < n; i++) {
        double tmp_angle = (double)i*2*PI/6.;
        positions.push_back({cos(tmp_angle), sin(tmp_angle), 0.});
      }
      positions[0] *= 1.5;
      positions[1] *= 5;
      positions[2] *= 10;
      positions[3] *= 5;
      positions[4] *= 1.5;
      positions[5] *= 1;
    }
    else if (std::strcmp(poly_str.c_str(), "rndbs 9-gon 1") == 0){
      n = 9;
      for (size_t i = 0; i < n-2; i++){
        greedy_triangulation.push_back({0, i+1, i+2});
      }
      for (int i = 0; i < n; i++) {
        double tmp_angle = (double)i*2*PI/(double)n;
        positions.push_back({cos(tmp_angle), sin(tmp_angle), 0.});
      }
      positions[0] *= 1;
      positions[1] *= 1.5;
      positions[2] *= 4;
      positions[3] *= 8;
      positions[4] *= 12;
      positions[5] *= 8;
      positions[6] *= 3;
      positions[7] *= 1.5;
      positions[8] *= 1;
    }
    else {
      throw std::runtime_error("no valid string provided\n");
    }
    // std::unique_ptr<ManifoldSurfaceMesh> poly_triangulated;
    // poly_triangulated.reset(new ManifoldSurfaceMesh(greedy_triangulation));
    mesh = new ManifoldSurfaceMesh(greedy_triangulation);
    geometry = new VertexPositionGeometry(*mesh);
    for (Vertex v : mesh->vertices()) {
        // Use the low-level indexers here since we're constructing
        printf("v %d\n", v.getIndex());
        geometry->inputVertexPositions[v] = positions[v.getIndex()];
    }
}

void update_solver(){
  forwardSolver = ForwardSolver(mesh, geometry, G);
  //assuming convex input here
  forwardSolver.hullMesh = new ManifoldSurfaceMesh(mesh->getFaceVertexList()); //mesh->copy().release();
  forwardSolver.hullGeometry = new VertexPositionGeometry(*forwardSolver.hullMesh); // geometry->copy().release();
  for (Vertex v: forwardSolver.hullMesh->vertices()){
    forwardSolver.hullGeometry->inputVertexPositions[v] = geometry->inputVertexPositions[v.getIndex()];
  }
  forwardSolver.compute_vertex_probabilities();


  // Register the mesh with polyscope
  psInputMesh = polyscope::registerSurfaceMesh(
      "input mesh",
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));
  draw_G();
}

void visualize_g_vec(){
  std::vector<Vector3> the_g_vec = {forwardSolver.current_g_vec};
  polyscope::PointCloudVectorQuantity *psG_vec = psG->addVectorQuantity("g_vec", the_g_vec);
  psG_vec->setEnabled(true);
  psG_vec->setVectorRadius(curve_radi_scale * 1.);
  psG_vec->setVectorLengthScale(0.2);
}

void visualize_contact_point(){
  // single point cloud for the single contact point
  if (forwardSolver.curr_state.first == forwardSolver.curr_state.second){
    std::vector<Vector3> curr_state_pos = {forwardSolver.hullGeometry->inputVertexPositions[forwardSolver.curr_state.first]}; // first and second should be the same since we just initialized.
    curr_state_pt = polyscope::registerPointCloud("current state", curr_state_pos);
    curr_state_pt->setEnabled(true);
    curr_state_pt->setPointRadius(pt_cloud_radi_scale/2.);
  }
  else {
      Vertex v1 = forwardSolver.curr_state.first, 
             v2 = forwardSolver.curr_state.second;
      std::vector<std::array<size_t, 2>> edgeInds;
      std::vector<Vector3> positions;
      edgeInds.push_back({0, 1});
      Vector3 p1 = forwardSolver.hullGeometry->inputVertexPositions[v1],
              p2 = forwardSolver.hullGeometry->inputVertexPositions[v2];
      positions.push_back(p1); positions.push_back(p2);
      auto curr_state_segment_net = polyscope::registerCurveNetwork("current contact edge", positions, edgeInds);
      curr_state_segment_net->setRadius(curve_radi_scale/1.5);
      curr_state_segment_net->setColor({0., 0., 1.});
      curr_state_segment_net->setEnabled(true);
  }
}

void initialize_state_vis(){
  visualize_g_vec();
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
        for (std::string tmp_str: all_polygon_items){ // This enables not having to have a const char* arr[]. Or maybe I'm just a noob.
            bool is_selected = (all_polygons_current_item == tmp_str.c_str()); // You can store your selection however you want, outside or inside your objects
            if (ImGui::Selectable(tmp_str.c_str(), is_selected)){ // selected smth
                all_polygons_current_item = tmp_str;
                generate_polygon_example(all_polygons_current_item);
                update_solver();
            }
            if (is_selected)
                ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
        }
        ImGui::EndCombo();
    }

  if (ImGui::Button("build boundary curves")) {
    visualize_boundary_curves();
  }
  if (ImGui::Button("visualize vertex probabilities")) {
    visualize_vertex_probabilities();
  }
  if (ImGui::SliderFloat("vertex radi scale", &pt_cloud_radi_scale, 0., 1.)) {
    draw_G();
  }
  if (ImGui::SliderFloat("edge radi scale", &curve_radi_scale, 0., 1.)) visualize_edge_probabilities();

  if (ImGui::Button("visualize edge probabilities")) {
    visualize_edge_probabilities();
  }

  if (ImGui::SliderFloat("G radi scale", &G_scale, -3., 3.)||
      ImGui::SliderFloat("G angle", &G_angle, 0., 2*PI)) {
    G = {cos(G_angle)*G_scale, sin(G_angle) * G_scale, 0};
    forwardSolver.G = G;
    draw_G();
    visualize_edge_probabilities();
  }

  if (ImGui::InputInt("compute empirical edge probabilities", &sample_count, 100)) {
    forwardSolver.build_next_edge_tracer();
    forwardSolver.empirically_build_probabilities(sample_count);
  }
  if (ImGui::SliderFloat("g-vec angle", &g_vec_angle, 0., 2*PI)) {
    initial_g_vec = {cos(g_vec_angle), sin(g_vec_angle), 0};
    forwardSolver.find_contact(initial_g_vec);
    forwardSolver.build_next_edge_tracer();
    initialize_state_vis();
  }
  if (ImGui::Button("reset state")) {
    forwardSolver.find_contact(initial_g_vec);
    forwardSolver.build_next_edge_tracer();
    initialize_state_vis();
  }
  if (ImGui::Button("next state")) {
    forwardSolver.next_state();
    visualize_g_vec();
    visualize_contact_point();
  }
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
  generate_polygon_example(all_polygons_current_item);
  G = {0.,0,0};
  update_solver();
  // build the solver
  
  // Load mesh
  // std::tie(mesh, geometry) = readManifoldSurfaceMesh(args::get(inputFilename));

  // Initialize polyscope
  polyscope::init();
  
  // Set the callback function
  polyscope::state::userCallback = myCallback;
  polyscope::view::style = polyscope::view::NavigateStyle::Planar;
  // Set vertex tangent spaces
  
  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
