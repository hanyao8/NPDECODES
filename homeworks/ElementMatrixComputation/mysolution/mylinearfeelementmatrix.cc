/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "mylinearfeelementmatrix.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>

namespace ElementMatrixComputation {

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix<double, 4, 4> MyLinearFEElementMatrix::Eval(
    const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  // Matrix storing corner coordinates in its columns
  auto vertices = geo_ptr->Global(ref_el.NodeCoords());
  // Matrix for returning element matrix
  Eigen::Matrix<double, 4, 4> elem_mat;

  //====================
  // Your code goes here
  //====================
  Eigen::Matrix<double,4,4> laplace_elem_mat;
  Eigen::Matrix<double,4,4> mass_elem_mat;

  lf::uscalfe::LinearFELaplaceElementMatrix laplace_elmat_builder;
  laplace_elem_mat = laplace_elmat_builder.Eval(cell);

  switch (ref_el) {
      case (lf::base::RefEl::kTria()): {
          double area = 
              0.5*
              ( (vertices(0,1)-vertices(0,0))*(vertices(1,2)-vertices(1,0))-
                (vertices(1,1)-vertices(1,0))*(vertices(0,2)-vertices(0,0)))
          mass_elem_mat <<
          2.0,1.0,1.0,0.0,
          1.0,2.0,1.0.0.0,
          1.0,1.0,2.0,0.0,
          0.0,0.0,0.0,0.0;

          mass_elem_mat *= area/12.0;
          break;
      }
      case (lf::base::RefEL::kQuad()): {
          double area = 
              (vertices(0,1)-vertices(0,0))*(vertices(1,3)-vertices(1,0));
          mass_elem_mat <<
          4.0,2.0,1.0,2.0,
          2.0,4.0,2.0,1.0,
          1.0,2.0,4.0,2.0,
          2.0,1.0,2.0,4.0;

          mass_elem_mat *= area/36.0;
          break;
      }
  }
  elem_mat = laplace_elem_mat + mass_elem_mat;
  return elem_mat;
}
/* SAM_LISTING_END_1 */
}  // namespace ElementMatrixComputation
