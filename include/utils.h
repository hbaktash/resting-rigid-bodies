#pragma once


#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <Eigen/Core>

// using namespace geometrycentral;
// using namespace geometrycentral::surface;


Eigen::Vector3d to_eigen(const geometrycentral::Vector3& _v);
geometrycentral::Vector3 to_geometrycentral(const Eigen::Vector3d& _v);
geometrycentral::Vector<double> tinyAD_flatten(geometrycentral::DenseMatrix<double> mat);
geometrycentral::DenseMatrix<double> unflat_tinyAD(geometrycentral::Vector<double> flat_mat);
geometrycentral::SparseMatrix<double> tinyADify_barrier_hess(std::vector<geometrycentral::DenseMatrix<double>> hessians);
geometrycentral::Vector<double> vec32vec(geometrycentral::Vector3 v);
geometrycentral::Vector3 vec_to_GC_vec3(geometrycentral::Vector<double> vec);
geometrycentral::DenseMatrix<double> vertex_data_to_matrix(geometrycentral::surface::VertexData<geometrycentral::Vector3> positions);
geometrycentral::surface::VertexData<geometrycentral::Vector3> vertex_matrix_to_data(Eigen::MatrixXd positions, geometrycentral::surface::ManifoldSurfaceMesh& mesh);
geometrycentral::DenseMatrix<double> face_data_to_matrix(geometrycentral::surface::FaceData<geometrycentral::Vector3> fdata);
// Eigen::SparseMatrix<double> tinyADify_constraint_mat(Eigen::MatrixXd A);
double binomial_dist(int n, int k);


// IO file
void write_mesh_obj_with_stability_material(geometrycentral::surface::SurfaceMesh& mesh, geometrycentral::surface::FaceData<double> prob,
                                            geometrycentral::surface::VertexPositionGeometry& geometry, std::string filename);




// // from GOAST


// // For each integer i in the mask, the i-th row-column of mat is set to e_i
// template<typename MatrixType>
// void applyMaskToMatrix( const std::vector<int>& mask, MatrixType& mat, bool setDiagonalOne = true  ) {

//   if ( setDiagonalOne )
//     for ( int mIdx : mask )
//       mat.coeffRef( mIdx, mIdx ) = 1.;

//   for ( int k = 0; k < mat.outerSize(); ++k ) {
//     for ( typename MatrixType::InnerIterator it( mat, k ); it; ++it ) {
//       if ( std::find_if( mask.begin(), mask.end(),
//                          [it]( const int &s ) { return (s == it.row() || s == it.col()); } ) != mask.end())
//         it.valueRef() = (it.row() == it.col() && setDiagonalOne ) ? 1. : 0.;
//     }
//   }
// }

// // For each integer i in the mask, the i-th entry of vec to 0
// template<typename VectorType>
// void applyMaskToVector( const std::vector<int>& mask, VectorType& vec ) {
//     // run over mask entries
//     for (uint k= 0; k<mask.size(); ++k){
//       int majorIdx = mask[k];      
//       // mask vector
//       if ( !(majorIdx < vec.size()) ){
//           std::cerr << "index = " << majorIdx << " whereas size = " << vec.size() << std::endl;
//           throw BasicException("applyMaskToVector() index out of bounds!");
//       }
//       vec[majorIdx] = 0.;
//     }
// }

// // For each integer i in the mask, the i-th row-column of mat is set to e_i, and the i-th entry of vec to 0
// template<typename MatrixType, typename VectorType>
// void applyMaskToSymmetricMatrixAndVector( const std::vector<int>& mask, MatrixType& mat, VectorType& vec ) {
//   applyMaskToVector<VectorType>( mask, vec );
//   applyMaskToSymmetricMatrix<MatrixType>(mask, mat );
// }

