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
#include <stan/math.hpp>
#include "forward3D.h"
#include "utils.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

// #include "visual_utils.h"

// Add missing include path for autodiff library
// #include <autodiff/reverse/var/eigen.hpp>
// #include <autodiff/reverse/var.hpp>

using namespace geometrycentral;
using namespace geometrycentral::surface;


std::vector<std::pair<Vector3, double>> 
normal_prob_assignment(std::string shape_name);

std::vector<std::pair<Vector3, double>> 
normal_prob_assignment_fair(Forward3DSolver *tmp_solver, size_t dice_side_count);

FaceData<double> 
manual_stable_only_face_prob_assignment(Forward3DSolver *tmp_solver, std::vector<std::pair<Vector3, double>> normal_prob_pairs);

std::vector<std::tuple<std::vector<Face>, double, Vector3>> 
manual_clustered_face_prob_assignment(Forward3DSolver *tmp_solver, std::vector<std::pair<Vector3, double>> normal_prob_pairs);

std::vector<std::pair<Vector3, double>> 
update_normal_prob_assignment(Forward3DSolver *tmp_solver,
                              std::vector<std::tuple<std::vector<Face>, double, Vector3>> clustered_face_normals,
                              bool take_max_prob_face);

// boundary of regions leading to different faces
class BoundaryNormal {
    public:
        static size_t counter;
        const size_t index;
        
        Vector3 normal;
        // autodiff::Vector3var normal_ad;

        std::vector<BoundaryNormal*> neighbors;

        Vertex host_v;
        Edge host_e;
        // "stable" face basins on attraction that this boundary separates; null for vertex equilibrias i.e. maxima
        Face f1, f2; // normals evaluated at area compute time only

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
        // ManifoldSurfaceMesh* mesh;
        // VertexPositionGeometry* geometry;
        
        // constructor
        BoundaryBuilder(Forward3DSolver *forward_solver_);

        // containers
        VertexData<BoundaryNormal*> vertex_boundary_normal;
        EdgeData<std::vector<BoundaryNormal*>> edge_boundary_normals; // could have multiple on a non-singular edge
        // EdgeData<std::vector<Vector3>> edge_boundary_normals;
        // FaceData<std::vector<BoundaryNormal*>> face_attraction_boundary; // for future passes over face region boundaries  
        FaceData<std::vector<std::tuple<BoundaryNormal*, BoundaryNormal*, double>>> face_chain_area;
        FaceData<double> face_region_area;

        // to store local minima distances to neighboring regions
        FaceData<FaceData<double>> MS_neighbor_distances;
        
        // FaceData<autodiff::Vector3var> face_normals_ad;

        FaceData<Eigen::MatrixX3d> df_dv_grads_ad;
        FaceData<Eigen::Vector3d> df_dG_grads;
        // FaceData<VertexData<Eigen::Vector3d>> df_dv_grads_ad;

        std::vector<Edge> find_terminal_edges();
        // backtrack and boundary normals starting from singular edges leading to different stable faces 
        void build_boundary_normals();
        // flow back from a edge with boundary normal; till u find a source
        bool flow_back_boundary_on_edge(BoundaryNormal* bnd_normal, Edge src_e, Vertex common_vertex,
                                        double f1_area_sign, double f2_area_sign);

        void print_area_of_boundary_loops();
        double get_fair_dice_energy(size_t side_count);

        double print_pairwise_distances(bool verbose = false);


        // TODO template
        template <typename Scalar>
        static Scalar dice_energy(Eigen::MatrixX3<Scalar> hull_positions, Eigen::Vector3<Scalar> G,
                                  Forward3DSolver &tmp_solver, double bary_reg, double co_planar_reg, double cluster_distance_reg, double unstable_attaction_thresh,
                                  std::string policy, std::vector<std::pair<Vector3, double>> normal_prob_assignment, 
                                  size_t side_count, bool verbose);
        template <typename Scalar>
        static Eigen::Vector3<Scalar> point_to_segment_normal(Eigen::Vector3<Scalar> P, Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B);
        template <typename Scalar>
        static Eigen::Vector3<Scalar> intersect_arcs(Eigen::Vector3<Scalar> v_normal, Eigen::Vector3<Scalar> R2, Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B);
        template <typename Scalar>
        static Scalar triangle_patch_signed_area_on_sphere(Eigen::Vector3<Scalar> A, Eigen::Vector3<Scalar> B, Eigen::Vector3<Scalar> C);

        // Regularizers
        //  --- for ManualCluster policy --- 
        template <typename Scalar>
        static Scalar 
        single_cluster_coplanar_e(std::vector<Face> faces, 
                                    FaceData<Eigen::Vector3<Scalar>> face_normals,
                                    FaceData<Face> face_last_face, double unstable_attaction_thresh, 
                                    bool verbose);
        
        template <typename Scalar>
        static Scalar 
        single_cluster_bary_e(std::vector<Face> faces, FaceData<Face> face_last_face,
                              Eigen::MatrixX3<Scalar> hull_positions, Eigen::Vector3<Scalar> G,
                              FaceData<Eigen::Vector3<Scalar>> face_normals,
                              FaceData<std::set<Vertex>> face_region_boundary_vertices);
        template <typename Scalar>
        static std::pair<Eigen::Vector3<Scalar>, Eigen::Vector3<Scalar>> 
        region_and_stables_barycenters(std::vector<Face> faces, FaceData<Face> face_last_face, 
                                       Eigen::MatrixX3<Scalar> hull_positions, Eigen::Vector3<Scalar> G,
                                       FaceData<Eigen::Vector3<Scalar>> face_normals,
                                       FaceData<std::set<Vertex>> face_region_boundary_vertices);

        // --- for Fair and Manual policies ---
        template <typename Scalar>
        static Scalar 
        single_DAG_cluster_bary_e(Face stable_face, Eigen::Vector3<Scalar> stable_face_normal,
                                Eigen::MatrixX3<Scalar> hull_positions, Eigen::Vector3<Scalar> G,
                                FaceData<std::set<Vertex>> face_region_boundary_vertices);
};

double hull_update_line_search(Eigen::MatrixX3d dfdv, Eigen::MatrixX3d hull_positions, Eigen::Vector3d G_vec, 
                               double bary_reg, double coplanar_reg, double cluster_distance_reg, double unstable_attaction_thresh,
                               std::string policy_general, std::vector<std::pair<Vector3, double>> normal_prob_assignment, 
                               size_t dice_side_count, double step_size, double decay, bool frozen_G, 
                               size_t max_iter, double step_tol);

#include "boundary_tools.impl.h"
