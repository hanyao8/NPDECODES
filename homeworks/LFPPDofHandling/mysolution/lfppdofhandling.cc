/**
 * @file
 * @brief NPDE homework "Handling degrees of freedom (DOFs) in LehrFEM++"
 * @author Julien Gacon
 * @date March 1st, 2019
 * @copyright Developed at ETH Zurich
 */

#include "lfppdofhandling.h"

#include <Eigen/Dense>
#include <array>
#include <memory>

#include "lf/assemble/assemble.h"
#include "lf/base/base.h"
#include "lf/geometry/geometry.h"
#include "lf/mesh/mesh.h"
#include "lf/mesh/utils/utils.h"

namespace LFPPDofHandling {

/* SAM_LISTING_BEGIN_1 */
std::array<std::size_t, 3> countEntityDofs(
    const lf::assemble::DofHandler &dofhandler) {
  std::array<std::size_t, 3> entityDofs;
  //====================
  // Your code goes here
  //====================
  for (int codim=0; codim<3; codim++){
    entityDofs[codim] = 0;
  }

  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();

  for (int codim=0; codim<3; codim++){
    for (const auto *entity : mesh->Entities(codim)){
      if (entity->RefEl() == lf::base::RefEl::kQuad()){
        throw("Error");
      }
      entityDofs[codim] += dofhandler.NumInteriorDofs(*entity);
    }
  }
  return entityDofs;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
std::size_t countBoundaryDofs(const lf::assemble::DofHandler &dofhandler) {
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  // given an entity, bd\_flags(entity) == true, if the entity is on the
  // boundary
  lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags(
      lf::mesh::utils::flagEntitiesOnBoundary(mesh));
  std::size_t no_dofs_on_bd = 0;
  //====================
  // Your code goes here
  //====================

  for (const auto *vertex : mesh->Entities(2)){
    if (bd_flags(*vertex)){
      no_dofs_on_bd+=1;
    }
  }
  return no_dofs_on_bd;
}
/* SAM_LISTING_END_2 */

// clang-format off
/* SAM_LISTING_BEGIN_3 */
double integrateLinearFEFunction(
    const lf::assemble::DofHandler& dofhandler,
    const Eigen::VectorXd& mu) {
  double I = 0;
  //====================
  // Your code goes here
  //====================

  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  double area;
  for (const auto *cell : mesh->Entities(0)){
    lf::base::size_type cell_local_dofs = dofhandler.NumLocalDofs(*cell);
    if (cell_local_dofs!=3){
      throw("Error");
    }
    lf::geometry::Geometry *cell_geo = cell->Geometry();
    area = lf::geometry::Volume(*cell_geo);
    const auto global_idxs = dofhandler.GlobalDofIndices(*cell);
    for (auto idx_p=global_idxs.begin() ; idx_p<global_idxs.end() ; ++idx_p){
      I += area/3.0 * mu(*idx_p);
    }
  }
  return I;
}
/* SAM_LISTING_END_3 */
// clang-format on

/* SAM_LISTING_BEGIN_4 */
double integrateQuadraticFEFunction(const lf::assemble::DofHandler &dofhandler,
                                    const Eigen::VectorXd &mu) {
  double I = 0;
  //====================
  // Your code goes here
  //====================
  //
  std::shared_ptr<const lf::mesh::Mesh> mesh = dofhandler.Mesh();
  double area;
  for (const auto *cell : mesh->Entities(0)){
    lf::base::size_type cell_local_dofs = dofhandler.NumLocalDofs(*cell);
    if (cell_local_dofs!=6){
      throw("Error");
    }
    lf::geometry::Geometry *cell_geo = cell->Geometry();
    area = lf::geometry::Volume(*cell_geo);
    const auto global_idxs = dofhandler.GlobalDofIndices(*cell);
    for (int i=3; i<6; i++){
      I += (area/3.0 * mu[global_idxs[i]]);
    }
  }
  return I;
}
/* SAM_LISTING_END_4 */

/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd convertDOFsLinearQuadratic(
    const lf::assemble::DofHandler &dofh_Linear_FE,
    const lf::assemble::DofHandler &dofh_Quadratic_FE,
    const Eigen::VectorXd &mu) {
  if (dofh_Linear_FE.Mesh() != dofh_Quadratic_FE.Mesh()) {
    throw "Underlying meshes must be the same for both DOF handlers!";
  }
  std::shared_ptr<const lf::mesh::Mesh> mesh =
      dofh_Linear_FE.Mesh();                          // get the mesh
  Eigen::VectorXd zeta(dofh_Quadratic_FE.NumDofs());  // initialise empty zeta
  // safety guard: always set zero if you're not sure to set every entry later
  // on for us this shouldn't be a problem, but just to be sure
  zeta.setZero();

  for (const auto *cell : mesh->Entities(0)) {
    // check if the spaces are actually linear and quadratic
    //====================
    // Your code goes here
    //====================
    // get the global dof indices of the linear and quadratic FE spaces, note
    // that the vectors obey the LehrFEM++ numbering, which we will make use of
    // lin\_dofs will have size 3 for the 3 dofs on the nodes and
    // quad\_dofs will have size 6, the first 3 entries being the nodes and
    // the last 3 the edges
    //====================
    // Your code goes here
    // assign the coefficients of mu to the correct entries of zeta, use
    // the previous subproblem 2-9.a
    //====================

    if (dofh_Linear_FE.NumLocalDofs(*cell) != 3 ||
        dofh_Quadratic_FE.NumLocalDofs(*cell) != 6){
      throw("Error");
    }
    nonstd::span<const lf::assemble::gdof_idx_t> lin_idxs = dofh_Linear_FE.GlobalDofIndices(*cell);
    nonstd::span<const lf::assemble::gdof_idx_t> quad_idxs = dofh_Quadratic_FE.GlobalDofIndices(*cell);
    for (int i=0; i<3; i++){
      zeta(quad_idxs[i]) = mu(lin_idxs[i]);
      zeta(quad_idxs[i+3]) = 0.5*mu(lin_idxs[i])+0.5*mu(lin_idxs[(i+1)%3]);
    }
  }
  return zeta;
}
/* SAM_LISTING_END_5 */

}  // namespace LFPPDofHandling
