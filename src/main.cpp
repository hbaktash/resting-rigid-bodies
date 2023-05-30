#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
// #include "polyscope/curve_network.h"
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


ForwardSolver forwardSolver;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psInputMesh, *dummy_psMesh1;

float pt_cloud_radi_scale = 1.;

void visualize_boundary_curves(){
  std::vector<std::vector<size_t>> dummy_face{{1,1,1}};
  dummy_psMesh1 = polyscope::registerSurfaceMesh(
      "dummy mesh1",
      geometry->inputVertexPositions, dummy_face);
  auto positions = forwardSolver.hullGeometry->inputVertexPositions.toVector();
  std::vector<std::array<size_t, 2>> edgeInds;
  for (Edge e: forwardSolver.hullMesh->edges()){
    if (e.isBoundary()){
      edgeInds.push_back({e.firstVertex().getIndex(), e.secondVertex().getIndex()});
    }
  }
  dummy_psMesh1->addSurfaceGraphQuantity("convex boundary", positions, edgeInds);
}


void visualize_vertex_probabilities(){
  
  for (Vertex v: forwardSolver.hullMesh->vertices()){
    std::vector<Vector3> positions = {forwardSolver.hullGeometry->inputVertexPositions[v]};
    polyscope::PointCloud* psCloud = polyscope::registerPointCloud("v" + std::to_string(v.getIndex()), positions);
    // set some options
    psCloud->setPointRadius(forwardSolver.vertex_probabilities[v] * pt_cloud_radi_scale);
    psCloud->setPointRenderMode(polyscope::PointRenderMode::Sphere);
  }

}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
  if (ImGui::Button("build boundary curves")) {
    visualize_boundary_curves();
  }
  if (ImGui::Button("visualize vertex probabilities")) {
    visualize_vertex_probabilities();
  }
  if (ImGui::SliderFloat("radi scale", &pt_cloud_radi_scale, 0., 1.)) visualize_vertex_probabilities();
  // ImGui::SliderFloat("param", &param1, 0., 100.);
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
    else {
      throw std::runtime_error("no valid string provided\n");
    }
    // std::unique_ptr<ManifoldSurfaceMesh> poly_triangulated;
    // poly_triangulated.reset(new ManifoldSurfaceMesh(greedy_triangulation));
    mesh = new ManifoldSurfaceMesh(greedy_triangulation);
    geometry = new VertexPositionGeometry(*mesh);
    for (Vertex v : mesh->vertices()) {
        // Use the low-level indexers here since we're constructing
        geometry->inputVertexPositions[v] = positions[v.getIndex()];
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
  generate_polygon_example("cube");
  // build the solver
  Vector3 G = {0,0,0};
  forwardSolver = ForwardSolver(mesh, geometry, G);
  //assuming convex input here
  forwardSolver.hullMesh = new ManifoldSurfaceMesh(mesh->getFaceVertexList()); //mesh->copy().release();
  forwardSolver.hullGeometry = new VertexPositionGeometry(*forwardSolver.hullMesh); // geometry->copy().release();
  for (Vertex v: forwardSolver.hullMesh->vertices()){
    forwardSolver.hullGeometry->inputVertexPositions[v] = geometry->inputVertexPositions[v.getIndex()];
  }
  forwardSolver.compute_vertex_probabilities();
  
  // Load mesh
  // std::tie(mesh, geometry) = readManifoldSurfaceMesh(args::get(inputFilename));

  // Initialize polyscope
  polyscope::init();
  // Register the mesh with polyscope
  psInputMesh = polyscope::registerSurfaceMesh(
      "input mesh",
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  
  // Set vertex tangent spaces
  
  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
