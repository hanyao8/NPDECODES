#include "incidencematrices.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <array>
#include <memory>

namespace IncidenceMatrices {

/** @brief Create the mesh consisting of a triangle and quadrilateral
 *         from the exercise sheet.
 * @return Shared pointer to the hybrid2d mesh.
 */
std::shared_ptr<lf::mesh::Mesh> createDemoMesh() {
  // builder for a hybrid mesh in a world of dimension 2
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);

  // Add points
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{0, 0});    // (0)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{1, 0});    // (1)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{1, 1});    // (2)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{0, 1});    // (3)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{0.5, 1});  // (4)

  // Add the triangle
  // First set the coordinates of its nodes:
  Eigen::MatrixXd nodesOfTria(2, 3);
  nodesOfTria << 1, 1, 0.5, 0, 1, 1;
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),  // we want a triangle
      std::array<lf::mesh::Mesh::size_type, 3>{
          {1, 2, 4}},  // indices of the nodes
      std::make_unique<lf::geometry::TriaO1>(nodesOfTria));  // node coords

  // Add the quadrilateral
  Eigen::MatrixXd nodesOfQuad(2, 4);
  nodesOfQuad << 0, 1, 0.5, 0, 0, 0, 1, 1;
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kQuad(),
      std::array<lf::mesh::Mesh::size_type, 4>{{0, 1, 4, 3}},
      std::make_unique<lf::geometry::QuadO1>(nodesOfQuad));

  std::shared_ptr<lf::mesh::Mesh> demoMesh_p = mesh_factory_ptr->Build();

  return demoMesh_p;
}

/** @brief Compute the edge-vertex incidence matrix G for a given mesh
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *        such as lf::mesh::hybrid2d::Mesh)
 * @return The edge-vertex incidence matrix as Eigen::SparseMatrix<int>
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<int> computeEdgeVertexIncidenceMatrix(
    const lf::mesh::Mesh &mesh) {
  // Store edge-vertex incidence matrix here
  Eigen::SparseMatrix<int, Eigen::RowMajor> G;
  std::cout<<"in compute G"<<std::endl;

  //====================
  // Your code goes here
  //====================

  const lf::mesh::Mesh::size_type numVertices = mesh.NumEntities(2);
  const lf::mesh::Mesh::size_type numEdges = mesh.NumEntities(1);

  G = Eigen::SparseMatrix<int, Eigen::RowMajor>(numEdges,numVertices);
  G.reserve(Eigen::VectorXi::Constant(numEdges,2));

  for (const lf::mesh::Entity *edge : mesh.Entities(1)){
    lf::mesh::Mesh::size_type edgeIdx = mesh.Index(*edge);
    auto vertices = edge->SubEntities(1);

    lf::mesh::Mesh::size_type firstvertexIdx = mesh.Index(*vertices[0]);
    lf::mesh::Mesh::size_type secondvertexIdx = mesh.Index(*vertices[1]);

    G.coeffRef(edgeIdx,firstvertexIdx) += 1;
    G.coeffRef(edgeIdx,secondvertexIdx) -= 1;
  }

  return G;
}
/* SAM_LISTING_END_1 */

/** @brief Compute the cell-edge incidence matrix D for a given mesh
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *        such as lf::mesh::hybrid2d::Mesh)
 * @return The cell-edge incidence matrix as Eigen::SparseMatrix<int>
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::SparseMatrix<int> computeCellEdgeIncidenceMatrix(
    const lf::mesh::Mesh &mesh) {
  // Store cell-edge incidence matrix here
  Eigen::SparseMatrix<int, Eigen::RowMajor> D;

  std::cout<<"in compute D"<<std::endl;

  //====================
  // Your code goes here
  //====================

  const lf::mesh::Mesh::size_type numEdges = mesh.NumEntities(1);
  const lf::mesh::Mesh::size_type numCells = mesh.NumEntities(0);

  D = Eigen::SparseMatrix<int, Eigen::RowMajor>(numCells,numEdges);
  D.reserve(Eigen::VectorXi::Constant(numCells,4));

  for (const lf::mesh::Entity *cell : mesh.Entities(0)){
    lf::mesh::Mesh::size_type cellIdx = mesh.Index(*cell);
    auto edges = cell->SubEntities(1);
    auto edgeOrientations = cell->RelativeOrientations(); //(Remember)

    auto edgeIt = edges.begin();
    auto orntIt = edgeOrientations.begin();
    for (;edgeIt!=edges.end() && orntIt!=edgeOrientations.end();
            ++edgeIt, ++orntIt) {
        lf::mesh::Mesh::size_type edgeIdx = mesh.Index(**edgeIt);
        D.coeffRef(cellIdx, edgeIdx) += lf::mesh::to_sign(*orntIt);
    }
  }
  return D;
}
/* SAM_LISTING_END_2 */

/** @brief For a given mesh test if the product of cell-edge and edge-vertex
 *        incidence matrix is zero: D*G == 0?
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *             such as lf::mesh::hybrid2d::Mesh)
 * @return true, if the product is zero and false otherwise
 */
/* SAM_LISTING_BEGIN_3 */
bool testZeroIncidenceMatrixProduct(const lf::mesh::Mesh &mesh) {
  bool isZero = false;

  //====================
  // Your code goes here
  //====================
  //

  Eigen::SparseMatrix<int> G = computeEdgeVertexIncidenceMatrix(mesh);
  Eigen::SparseMatrix<int> D = computeCellEdgeIncidenceMatrix(mesh);

  Eigen::MatrixXi Gd = Eigen::MatrixXi(G);
  Eigen::MatrixXi Dd = Eigen::MatrixXi(D);
  int a,b,c,d;
  a = Gd.rows();
  b = Gd.cols();
  c = Dd.rows();
  d = Dd.cols();
  std::cout<<a<<b<<c<<d<<std::endl;
  std::cout<<Gd<<std::endl;
  std::cout<<Dd<<std::endl;

  Eigen::SparseMatrix<int> DG = D*G;
  if (DG.norm() == 0){
    isZero = true;
  }
  return isZero;
}
/* SAM_LISTING_END_3 */

}  // namespace IncidenceMatrices
