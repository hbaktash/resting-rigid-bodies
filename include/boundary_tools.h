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
        VertexData<BoundaryNormal*> vertex_boundary_normals;
        EdgeData<BoundaryNormal*> edge_boundary_normals;

        // backtrack and boundary normals starting from singular edges leading to different stable faces 
        void build_boundary_normals();

        
};
