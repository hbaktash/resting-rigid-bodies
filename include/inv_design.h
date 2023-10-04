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
#include "boundary_tools.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


class InverseSolver{
    public:
        Forward3DSolver* forwardSolver;
        BoundaryBuilder* boundaryBuilder;

        FaceData<double> goal_area;

        InverseSolver(){}
        InverseSolver(BoundaryBuilder* boundaryBuilder);

        // distribution goals
        void set_fair_distribution();

        // gradient computation
        FaceData<Vector3> per_face_G_gradient;
        void find_per_face_G_grads();
        Vector3 find_total_g_grad();
};