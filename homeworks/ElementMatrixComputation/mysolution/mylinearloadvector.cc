/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler, edited by Oliver Rietmann
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "mylinearloadvector.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <functional>

namespace ElementMatrixComputation {

namespace {

/* SAM_LISTING_BEGIN_1 */
Eigen::Vector4d computeLoadVector(
    const Eigen::MatrixXd &vertices,
    std::function<double(const Eigen::Vector2d &)> f) {
  // Number of nodes of the element: triangles = 3, rectangles = 4
  const int num_nodes = vertices.cols();
  // Vector for returning element vector
  Eigen::Vector4d elem_vec = Eigen::Vector4d::Zero();

  //====================
  // Your code goes here
  //====================
 
  Eigen::MatrixXd midpoints(2,num_nodes);
  
  double area;
  switch(num_nodes){
    case 3: {
      area = 0.5*(vertices(0,1)-vertices(0,0)*(vertices(1,2)-vertices(1,0))-
                  vertices(0,2)-vertices(0,0)*(vertices(1,1)-vertices(1,0))
              );
      midpoints<<
        vertices(0,0)+vertices(0,1),
        vertices(0,1)+vertices(0,2),
        vertices(0,2)+vertices(0,0),
        vertices(1,0)+vertices(1,1),
        vertices(1,1)+vertices(1,2),
        vertices(1,2)+vertices(1,0);
      break;
    }
    case 4: {
      area = (vertices(0,1)-vertices(0,0))*(vertices(1,3)-vertices(1,0));
      midpoints<<
        vertices(0,0)+vertices(0,1),
        vertices(0,1)+vertices(0,2),
        vertices(0,2)+vertices(0,3),
        vertices(0,3)+vertices(0,0),
        vertices(1,0)+vertices(1,1),
        vertices(1,1)+vertices(1,2),
        vertices(1,2)+vertices(1,3),
        vertices(1,3)+vertices(1,0);
      break;
    }
  }
  midpoints *= 0.5;

  Eigen::VectorXd fvals = Eigen::VectorXd::Zero(num_nodes);
  for (int i=0;i<num_nodes;i++){
    fvals[i] = f(midpoints.col(i));
  }

  for (int i=0;i<num_nodes;i++){
    elem_vec[i%num_nodes] += 0.5*fvals[i];
    elem_vec[(i+1)%num_nodes] += 0.5*fvals[i];
  }
  elem_vec *= area/(double)num_nodes;
  return elem_vec;
}
/* SAM_LISTING_END_1 */

}  // namespace

Eigen::Vector4d MyLinearLoadVector::Eval(const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  const lf::base::size_type num_nodes{ref_el.NumNodes()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();

  // Matrix storing corner coordinates in its columns
  auto vertices = geo_ptr->Global(ref_el.NodeCoords());

  return computeLoadVector(vertices, f_);
}

}  // namespace ElementMatrixComputation
