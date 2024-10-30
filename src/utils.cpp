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

#include "utils.h"



// vector stuff
Eigen::Vector3d to_eigen(const geometrycentral::Vector3& _v){
    return Eigen::Vector3d(_v.x, _v.y, _v.z);
}

geometrycentral::Vector3 to_geometrycentral(
        const Eigen::Vector3d& _v) {
    return geometrycentral::Vector3 { _v.x(), _v.y(), _v.z() };
}

geometrycentral::Vector<double> tinyAD_flatten(geometrycentral::DenseMatrix<double> mat){
    size_t n = mat.rows();
    assert(mat.cols() == 3);
    geometrycentral::Vector<double> ans(n*3);
    for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j < 3; j++)
            ans(3*i + j) = mat(i,j);
    }
    return ans;
}

geometrycentral::DenseMatrix<double> unflat_tinyAD(geometrycentral::Vector<double> flat_mat){
    size_t n = flat_mat.size()/3;
    geometrycentral::DenseMatrix<double> mat(n, 3);
    for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j < 3; j++)
            mat.coeffRef(i,j) = flat_mat(3*i + j);
    }
    return mat;
}

geometrycentral::SparseMatrix<double> tinyADify_barrier_hess(std::vector<geometrycentral::DenseMatrix<double>> hessians){
    size_t n = hessians.size();
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(9*n);

    for (size_t i = 0; i < n; i++){
        geometrycentral::DenseMatrix<double> hess_i = hessians[i];
        for (size_t k = 0; k < 3; k++){
            for (size_t l = 0; l < 3; l++)
                tripletList.emplace_back(3 * i + k, 3 * i + l, hess_i(k, l));
        }
    }
    Eigen::SparseMatrix<double> tinyADfied_hess(3*n, 3*n);
    tinyADfied_hess.setFromTriplets(tripletList.begin(), tripletList.end());
    return tinyADfied_hess;
}


geometrycentral::Vector<double> vec32vec(geometrycentral::Vector3 v){
    geometrycentral::Vector<double> ans(3);
    ans[0] = v.x;
    ans[1] = v.y;
    ans[2] = v.z;
    return ans;
}

geometrycentral::Vector3 vec_to_GC_vec3(geometrycentral::Vector<double> vec){
    geometrycentral::Vector3 ans;
    ans.x = vec[0];
    ans.y = vec[1];
    ans.z = vec[2];
    return ans;
}

geometrycentral::DenseMatrix<double> vertex_data_to_matrix(geometrycentral::surface::VertexData<geometrycentral::Vector3> positions){
    size_t n = positions.getMesh()->nVertices();
    geometrycentral::DenseMatrix<double> mat(n, 3);
    for (geometrycentral::surface::Vertex v: positions.getMesh()->vertices()){
        geometrycentral::Vector3 p = positions[v];
        mat.row(v.getIndex()) = vec32vec(p);
    }
    return mat;
}


geometrycentral::DenseMatrix<double> face_data_to_matrix(geometrycentral::surface::FaceData<geometrycentral::Vector3> fdata){
    size_t n = fdata.getMesh()->nFaces();
    geometrycentral::DenseMatrix<double> mat(n, 3);
    for (geometrycentral::surface::Face f: fdata.getMesh()->faces()){
        geometrycentral::Vector3 p = fdata[f];
        mat.row(f.getIndex()) = vec32vec(p);
    }
    return mat;
}


geometrycentral::surface::VertexData<geometrycentral::Vector3> vertex_matrix_to_data(Eigen::MatrixXd positions, 
                                                                                     geometrycentral::surface::ManifoldSurfaceMesh& mesh){
    geometrycentral::surface::VertexData<geometrycentral::Vector3> ans(mesh);
    for (geometrycentral::surface::Vertex v: mesh.vertices()){
        ans[v] = vec_to_GC_vec3(positions.row(v.getIndex()));
    }
    return ans;
}

double binomial_dist(int n, int k) {
    if (k > n - k)  // Take advantage of symmetry
        k = n - k;
    double result = 1.0;
    for (int i = 1; i <= k; ++i) {
        result *= (double)(n - k + i)/(double)(2.*i);
    }
    result *= std::ldexp(1.0, -(n - k));
    return result;
}