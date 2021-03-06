/**
 * @ file boundarylength.cc
 * @ brief NPDE homework LengthOfBoundary code
 * @ author Christian Mitsch
 * @ date 03.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "boundarylength.h"

#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>

namespace LengthOfBoundary {

/* SAM_LISTING_BEGIN_1 */
double volumeOfDomain(const std::shared_ptr<lf::mesh::Mesh> mesh_p) {
  double volume = 0.0;
  //====================
  // Your code goes here
  //====================

//  for (const lf::mesh::Entity *cell : mesh_p->Entities){
//
  for (const lf::mesh::Entity *cell : mesh_p->Entities(0)){
    lf::geometry::Geometry *cell_geo = cell->Geometry();
    volume += lf::geometry::Volume(*cell_geo);
  }
  
  return volume;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double lengthOfBoundary(const std::shared_ptr<lf::mesh::Mesh> mesh_p) {
  double length = 0.0;
  //====================
  // Your code goes here
  //====================

  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p,1)};
  //returns a mapping 
  //that takes in *cell ptr object
  //and the mapping returns a boolean
  for (const lf::mesh::Entity *edge : mesh_p->Entities(1)){
    if (bd_flags(*edge)){
      lf::geometry::Geometry *edge_geo = edge->Geometry();
      length += lf::geometry::Volume(*edge_geo);
    }
  }

  return length;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
std::pair<double, double> measureDomain(std::string filename) {
  double volume, length;

  //====================
  // Your code goes here
  //====================

  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory),filename);
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader.mesh();

  volume = volumeOfDomain(mesh_p);
  length = lengthOfBoundary(mesh_p);

  return {volume, length};
}
/* SAM_LISTING_END_3 */

}  // namespace LengthOfBoundary
