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
// #include "geometrycentral/utilities/eigen_interop_helpers.h"
#include "boundary_tools.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "deformation.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;




class InverseSolver{
    public:
        Forward3DSolver* forwardSolver;
        BoundaryBuilder* boundaryBuilder;

        FaceData<double> goal_prob;

        InverseSolver(){}
        InverseSolver(BoundaryBuilder* boundaryBuilder);

        // distribution goals
        void set_fair_distribution();
        void set_fair_distribution_for_sink_faces();
        // smarter distribution update
        std::vector<Vector3> old_stable_normals;
        void update_fair_distribution(double normal_threshold);
        
        // // for keeping track of interior points
        // std::vector<Vector3> old_normals; 
        // std::vector<std::pair<Vector3, double>> old_to_new_normal_rotations;
        // VertexData<Face> interior_v_to_hull_f;
        // VertexData<double> interior_v_to_hull_f_hit_ratio;

        // void initialize_interior_vertex_trackers();
        // gradient computation; assuming regularity
        // note: pf = face region area
        // G grad; vertices frozen
        FaceData<Vector3> d_pf_d_G;            // convex-hull related
        void find_d_pf_d_Gs(bool check_FD = false);
        // accumulate over all faces
        Vector3 find_total_g_grad();
        // vertex grad; G frozen
        FaceData<VertexData<Vector3>> d_pf_dv; // convex-hull related
        void find_d_pf_dvs(bool check_FD = false); 
        // accumulate over all faces
        VertexData<Vector3> find_total_vertex_grads();

        // Uniform mass; G is dependent of Geometry
        // DG/dv
        VertexData<DenseMatrix<double>> dG_dv; // interior vertex related
        void find_dG_dvs();                    
        // dp_f/dv
        bool compute_global_G_effect = true;
        FaceData<VertexData<Vector3>> uni_mass_d_pf_dv;
        // after chain and product rule
        void find_uni_mass_d_pf_dv(bool check_FD = false);
        // accumulate over all faces
        FaceData<Face> flow_structure;
        VertexData<Vector3> find_uni_mass_total_vertex_grads(bool with_flow_structure = false,
                                                             double stable_normal_update_thres = -1.);
        
        // pre-step for updates 
        void subdivide_for_aggressive_updates(VertexData<Vector3> hull_updates);
        
        // position update for interior vertices
        // ARAP
        size_t arap_max_iter = 10;
        PositiveDefiniteSolver<double> *constrained_L_solver;
        geometrycentral::SparseMatrix<double> initial_Ls;
        geometrycentral::SparseMatrix<double> current_Ls;
        VertexPositionGeometry* initial_geometry;
        void save_initial_pos_and_Ls();
        bool use_old_Ls = false, 
             use_old_geometry = false, 
             libigl_ARAP = true;
        void update_Ls();
        DenseMatrix<double> solve_constrained_Laplace(Vector<size_t> interior_indices, DenseMatrix<double> old_pos, DenseMatrix<double> new_pos,
                                                      bool update_solver_decomp = true);
        VertexData<DenseMatrix<double>> find_rotations(DenseMatrix<double> old_pos, DenseMatrix<double> new_pos);
        VertexData<Vector3> ARAP_update_positions(VertexData<Vector3> hull_updates);
        VertexData<Vector3> laplace_update_positions(VertexData<Vector3> hull_updates);
        // Laplacian edit
        // Greedy stuff
        VertexData<Vector3> trivial_update_positions(VertexData<Vector3> hull_updates);
        VertexData<Vector3> diffusive_update_positions(VertexData<Vector3> hull_updates);

};