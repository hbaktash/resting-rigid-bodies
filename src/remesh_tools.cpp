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
    
    // Edge split callback
    auto funcPreSplit = [&](Edge oldE, double tSplit) -> Vector3 { 
      Vertex oldVA = oldE.firstVertex();
      Vertex oldVB = oldE.secondVertex();
      return (1. - tSplit) * ref_geometry->inputVertexPositions[oldVA] + tSplit * ref_geometry->inputVertexPositions[oldVB];
    };
    auto funcPostSplit = [&](Halfedge newHe1, Halfedge newHe2, double tSplit, Vector3 new_pos) -> void {
        ref_geometry->inputVertexPositions[newHe1.vertex()] = new_pos;
    };
    deformed_mm.registerEdgeSplitHandlers(funcPreSplit, funcPostSplit);
    
    // Edge collapse callback
    auto funcPreCollapse = [&](Edge oldE, double tSplit) -> Vector3 {
      Vertex oldVA = oldE.firstVertex();
      Vertex oldVB = oldE.secondVertex();
      return (1. - tSplit) * ref_geometry->inputVertexPositions[oldVA] + tSplit * ref_geometry->inputVertexPositions[oldVB];
    };
    auto funcPostCollapse = [&](Vertex new_v, double tSplit, Vector3 new_pos) -> void {
        ref_geometry->inputVertexPositions[new_v] = new_pos;
    };
    deformed_mm.registerEdgeCollapseHandlers(funcPreCollapse, funcPostCollapse);
    
    remesh(*mesh, *deformed_geometry, deformed_mm, options);
    
    mesh->compress();
    
    // debug visuals
    Vector3 vis_shift = {2,0,0};
    auto old_geo_psmesh = polyscope::registerSurfaceMesh("post remesh old geo", ref_geometry->inputVertexPositions + vis_shift, mesh->getFaceVertexList());
    auto new_geo_psmesh = polyscope::registerSurfaceMesh("post remesh def geo", deformed_geometry->inputVertexPositions, mesh->getFaceVertexList());
    old_geo_psmesh->setSurfaceColor({0.2,0.2,0.7});
    new_geo_psmesh->setSurfaceColor({0.7,0.2,0.2});
    polyscope::show();
}