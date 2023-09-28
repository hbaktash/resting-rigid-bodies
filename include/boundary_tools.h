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

#pragma once

// #include "geometrycentral/surface/manifold_surface_mesh.h"
// #include "geometrycentral/surface/vertex_position_geometry.h"
// #include "geometrycentral/surface/surface_point.h"
#include "forward3D.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


// boundary of regions leading to different faces
class BoundaryNormal {
    public:
        static size_t counter;
        const size_t index;
        
        Vector3 normal;
        std::vector<BoundaryNormal*> neighbors;

        Vertex host_v;
        Edge host_e;
        // faces that this boundary contributes to
        Face f1, f2; // null for maxima's (vertices)

        // constructors
        BoundaryNormal(): index(counter++){}
        BoundaryNormal(Vector3 _normal);

        // 
        void add_neighbor(BoundaryNormal* _neigh);
};


class BoundaryBuilder {
    public:
        Forward3DSolver* forward_solver;
        // assuming convexity
        // "same" pointer as the one in forward solver; here for easier access
        ManifoldSurfaceMesh* mesh;
        VertexPositionGeometry* geometry;
        
        // constructor
        BoundaryBuilder(Forward3DSolver *forward_solver_);

        // containers
        VertexData<BoundaryNormal*> vertex_boundary_normal;
        EdgeData<std::vector<BoundaryNormal*>> edge_boundary_normals; // could have multiple on a non-singular edge
        FaceData<std::vector<BoundaryNormal*>> face_attraction_boundary;  
        FaceData<std::vector<std::tuple<BoundaryNormal*, BoundaryNormal*, double>>> face_chain_area;
        FaceData<double> face_region_area;

        // backtrack and boundary normals starting from singular edges leading to different stable faces 
        void build_boundary_normals();

        // flow back from a edge with boundary normal; till u find a source
        void flow_back_boundary_on_edge(BoundaryNormal* bnd_normal, Edge src_e, Vertex common_vertex,
                                        double f1_area_sign, Vector3 f1_normal, Vector3 f2_normal);

        void print_area_of_boundary_loops();
        
};
