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

#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"

#include <remesh_tools.h>



void joint_remesh(ManifoldSurfaceMesh* mesh, VertexPositionGeometry *ref_geometry, VertexPositionGeometry *deformed_geometry, double target_edge_len){
    // not doint "joint" atm 
    RemeshOptions options = defaultRemeshOptions;
    options.targetEdgeLength = target_edge_len;
    
    MutationManager deformed_mm(*mesh, *deformed_geometry);
    // DEBUG
    Vector3 vis_shift = {2,0,0};
    // Edge split callback
    auto funcPreSplit = [&](Edge oldE, double tSplit) -> Vector3 { 
      Vertex oldVA = oldE.firstVertex();
      Vertex oldVB = oldE.secondVertex();
      // printf("pre spliting t: %f\n", tSplit);
      return (1. - tSplit) * ref_geometry->inputVertexPositions[oldVA] + tSplit * ref_geometry->inputVertexPositions[oldVB];
    };
    auto funcPostSplit = [&](Halfedge newHe1, Halfedge newHe2, double tSplit, Vector3 new_pos) -> void {
        ref_geometry->inputVertexPositions[newHe1.vertex()] = new_pos;
    };
    deformed_mm.registerEdgeSplitHandlers(funcPreSplit, funcPostSplit);
    
    // Edge collapse callback
    auto funcPreCollapse = [&](Edge oldE, double tSplit) -> Vector3 {
      Vertex oldVA = oldE.firstVertex(),
             oldVB = oldE.secondVertex();
      Vector3 pA = ref_geometry->inputVertexPositions[oldVA],
              pB = ref_geometry->inputVertexPositions[oldVB];
      // printf("pre collapsing t: %f\n", tSplit);
      tSplit = tSplit < 0. ? 0.5 : tSplit;
      Vector3 pMid = (1. - tSplit) * pA + tSplit * pB;
      // auto old_geo_psmesh = polyscope::registerSurfaceMesh("post remesh old geo", ref_geometry->inputVertexPositions + vis_shift, mesh->getFaceVertexList());
      // old_geo_psmesh->setSurfaceColor({0.3,0.3,0.8});
      // auto tmp_pc = polyscope::registerPointCloud("old edge", std::vector<Vector3>({pA + vis_shift, pB + vis_shift}));
      // tmp_pc->setPointColor({0.1, 0.8, 0.1});
      // auto tmp_pc2 = polyscope::registerPointCloud("old edge midP", std::vector<Vector3>({pMid + vis_shift}));
      // tmp_pc2->setPointColor({0.8, 0.1, 0.1});
      // polyscope::show();
      return pMid;
    };
    auto funcPostCollapse = [&](Vertex new_v, double tSplit, Vector3 new_pos) -> void {
        ref_geometry->inputVertexPositions[new_v] = new_pos;
    };
    deformed_mm.registerEdgeCollapseHandlers(funcPreCollapse, funcPostCollapse);
    
    EdgeData<bool> colapsibleEdges(*mesh),
                   flipableEdges(*mesh);
    VertexData<bool> repositionableVertices(*mesh);
    deformed_mm.setCollapsibleEdges(colapsibleEdges, false);
    deformed_mm.setFlippableEdges(flipableEdges, false);
    deformed_mm.setRepositionableVertices(repositionableVertices, false);
    remesh(*mesh, *deformed_geometry, deformed_mm, options);
    
    mesh->compress();
    
    // debug visuals
    auto old_geo_psmesh = polyscope::registerSurfaceMesh("post remesh old geo", ref_geometry->inputVertexPositions + vis_shift, mesh->getFaceVertexList());
    // auto new_geo_psmesh = polyscope::registerSurfaceMesh("post remesh def geo", deformed_geometry->inputVertexPositions, mesh->getFaceVertexList());
    old_geo_psmesh->setSurfaceColor({0.2,0.2,0.7});
    old_geo_psmesh->setEdgeWidth(1.);
    // new_geo_psmesh->setSurfaceColor({0.7,0.2,0.2});
    // polyscope::show();
}


bool split_only_remesh(ManifoldSurfaceMesh* mesh, VertexPositionGeometry *ref_geometry, VertexPositionGeometry *deformed_geometry, double target_edge_len){
  Vector3 vis_shift = {2,0,0};

  bool change_flag = false; 
  std::vector<Edge> toSplit;
  for (Edge e : mesh->edges()) {
    toSplit.push_back(e);
  }

  // actually splitting
  while (!toSplit.empty()) {
    Edge e = toSplit.back();
    toSplit.pop_back();
    double length_e = deformed_geometry->edgeLength(e);
    double threshold = 2. * target_edge_len;
    if (length_e > threshold) {
      Vector3 def_newPos = 0.5 * (deformed_geometry->inputVertexPositions[e.firstVertex()] + deformed_geometry->inputVertexPositions[e.secondVertex()]),
              ref_newPos = 0.5 * (ref_geometry->inputVertexPositions[e.firstVertex()] + ref_geometry->inputVertexPositions[e.secondVertex()]);
      Halfedge he = mesh->splitEdgeTriangular(e);
      change_flag = true;
      deformed_geometry->inputVertexPositions[he.vertex()] = def_newPos;
      ref_geometry->inputVertexPositions[he.vertex()] = ref_newPos;
    }
  }
  mesh->compress();
  // debug visuals
  auto old_geo_psmesh = polyscope::registerSurfaceMesh("post remesh old geo", ref_geometry->inputVertexPositions + vis_shift, mesh->getFaceVertexList());
  // auto new_geo_psmesh = polyscope::registerSurfaceMesh("post remesh def geo", deformed_geometry->inputVertexPositions, mesh->getFaceVertexList());
  old_geo_psmesh->setSurfaceColor({0.2,0.2,0.7});
  old_geo_psmesh->setEdgeWidth(1.);
  // new_geo_psmesh->setS
  return change_flag;
}