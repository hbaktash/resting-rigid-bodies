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
 #include "deformation.h"


geometrycentral::DenseMatrix<double> get_ARAP_positions(
                                       geometrycentral::DenseMatrix<double> old_pos_mat,
                                       geometrycentral::DenseMatrix<double> new_pos_mat, 
                                       geometrycentral::DenseMatrix<double> init_sol, 
                                       ManifoldSurfaceMesh &inner_mesh,
                                       geometrycentral::Vector<int> hull_indices){
    Eigen::MatrixXd V,U, bc;
    Eigen::MatrixXi F;
    Eigen::VectorXi S,b;
    igl::ARAPData arap_data;
    arap_data.max_iter = 20;
    SparseMatrix<double> L;
    F = inner_mesh.getFaceVertexMatrix<size_t>().cast<int>();
    V = old_pos_mat;
    b = hull_indices;
    // igl::arap_precomputation(old_pos_mat, F, 3, hull_indices, arap_data);
    igl::arap_precomputation(V, F, 3, b, arap_data);
    // geometrycentral::DenseMatrix<double> U = init_sol; 
    U = init_sol;
    bc = new_pos_mat(hull_indices, Eigen::all);
    igl::arap_solve(bc, arap_data, U);
    return U;
}