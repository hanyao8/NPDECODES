/**
 * @file
 * @brief NPDE homework ProjectionOntoGradients code
 * @author ?, Philippe Peter
 * @date December 2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <utility>

namespace ProjectionOntoGradients {

/* SAM_LISTING_BEGIN_1 */
class ElementMatrixProvider {
 public:
  Eigen::Matrix3d Eval(const lf::mesh::Entity &entity);
  bool isActive(const lf::mesh::Entity & /*entity*/) const { return true; }
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::Matrix3d ElementMatrixProvider::Eval(const lf::mesh::Entity &entity) {
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  const lf::geometry::Geometry *geo_ptr = entity.Geometry();
  Eigen::Matrix3d loc_mat;

  //====================
  // Your code goes here
  //====================
  double area;
  area = lf::Geometry::Volume(*geo_ptr);
  
  Eigen::Matrix3d lin_params;
  Eigen::Matrix<double,2,3> cell_corners = lf::Geometry::Corners(*geo_ptr);
  Eigen::Matrix3d lin_data;
  lin_data.col(0) = Eigen::VectorXd::Ones(3);
  lin_data.block<3,2>(0,1) = cell_corners.transpose();
  lin_params = (lin.data.inverse()).block<2,3>(1,0);

  loc_mat = area*(lin_params.transpose)*lin_params;

  return loc_mat;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR>
class GradProjRhsProvider {
 public:
  explicit GradProjRhsProvider(FUNCTOR f) : f_(f) {}

  Eigen::Vector3d Eval(const lf::mesh::Entity &entity);
  bool isActive(const lf::mesh::Entity & /*entity*/) const { return true; }

 private:
  FUNCTOR f_;
};
/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
template <typename FUNCTOR>
Eigen::Vector3d GradProjRhsProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity &entity) {
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  const lf::geometry::Geometry *geo_ptr = entity.Geometry();
  Eigen::Vector3d loc_vec;

  //====================
  // Your code goes here
  //====================

  Eigen::Matrix<double,2,3> G;
  Eigen::Matrix<double,2,3> cell_corners;

  cell_corners = lf::geometry::Corners(*geo_ptr);
  Eigen::Matrix3d lin_data;
  lin_data.col(0) = Eigen::VectorXd::Ones(3);
  lin_data.block<3,2>(0,1) = cell_corners.transpose();
  G = (lin_data.inverse()).block<2,3>(1,0);

  double area;
  area = lf::geometry::Volume(*geo_ptr);

  Eigen::Vector2d midpoint = (cell_corners.col(0)+
      cell_corners.col(1)+cell_corners.col(2))/3.0;
  Eigen::Vector2d f_mid = f_(midpoint);

  loc_vec = area*(G.transpose())*f_mid;

  return loc_vec;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
template <typename FUNCTOR>
Eigen::VectorXd projectOntoGradients(const lf::assemble::DofHandler &dofh,
                                     FUNCTOR f) {
  const lf::assemble::size_type N_dofs = dofh.NumDofs();
  Eigen::VectorXd sol_vec;

  // I. Build the (full) Galerkin matrix for the linear system.
  //====================
  // Your code goes here
  //====================

  std::shared_ptr<const lf::mesh::Mesh> mesh = dofh.Mesh();
  //cell oriented assembly
  lf::assemble::COOMatrix<double> A_coo(N_dofs,N_dofs);
  ProjectionOntoGradients::ElementMatrixProvider locmat_provider;
  lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
          0,dofh,dofh,locmat_provider,A_coo);

  // II. Build the (full) right hand side vector
  //====================
  // Your code goes here
  //====================
  Eigen::VectorXd phi = Eigen::VectorXd::Zero(N_dofs);
  ProjectionOntoGradients::GradProjRhsProvider locvec_provider;
  lf::assemble::AssembleVectorLocally(0,dofh,locvec_provider);

  // III. Enforce homogeneous dirichlet boundary conditions
  //====================
  // Your code goes here
  //====================

  //lf::mesh::utils::CodimMeshDataSet
  

  // IV. Solve the LSE using an Eigen solver
  //====================
  // Your code goes here
  //====================
  return sol_vec;
}
/* SAM_LISTING_END_5 */

}  // namespace ProjectionOntoGradients
