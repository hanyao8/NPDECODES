/**
 * @ file LinearFE1D.h
 * @ brief NPDE homework LinearFE1D code
 * @ author Christian Mitsch, Amélie Loher, Erick Schulz
 * @ date 11.11.2019
 * @ copyright Developed at ETH Zurich
 */

#include <Eigen/SparseCholesky>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

namespace LinearFE1D {

// Calculate the matrix entries corresponding to the Galerkin matrix
// of the Laplace operator with non-constant coefficient function alpha
// using composite midpoint integration rule. The entries are stored in
// the Eigen triplet data structure.
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR1>
std::vector<Eigen::Triplet<double>> computeA(const Eigen::VectorXd &mesh,
                                             FUNCTOR1 &&alpha) {
  // Nodes are indexed as 0=x_0 < x_1 < ... < x_N = 1
  unsigned N = mesh.size() - 1;
  // Initializing the vector of triplets whose size corresponds to the
  // number of entries in a (N+1) x (N+1) tridiagonal (band) matrix
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(3 * N + 1);

  //====================
  // Your code goes here
  //====================

  return triplets;
}  // computeA
/* SAM_LISTING_END_1 */

// Calculate the matrix entries corresponding to the mass matrix
// associated to the L2 pairing weighted by the variable coefficient
// function gamma using the trapezoidal integration rule.
/* SAM_LISTING_BEGIN_2 */
template <typename FUNCTOR1>
std::vector<Eigen::Triplet<double>> computeM(const Eigen::VectorXd &mesh,
                                             FUNCTOR1 &&gamma) {
  // Nodes are indexed as 0=x_0 < x_1 < ... < x_N = 1
  unsigned N = mesh.size() - 1;

  // The basis tent functions spanning the space of continuous piecewise
  // linear polynomial satisfy the cardinal basis property with respect
  // to the nodes of the mesh. Using  the trapezoidal rule to approximate
  // the integrals therefore lead to a diagonal (N+1) x (N+1) matrix
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(N + 1);

  //====================
  // Your code goes here
  //====================

  return triplets;
}  // computeM
/* SAM_LISTING_END_2 */

// Calculate the entries of the right hand side vector using
// the composite trapezoidal integration rule
/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR1>
Eigen::VectorXd computeRHS(const Eigen::VectorXd &mesh, FUNCTOR1 &&f) {
  // Nodes are indexed as 0=x_0 < x_1 < ... < x_N = 1
  unsigned N = mesh.size() - 1;
  // Initializing right hand side vector
  Eigen::VectorXd rhs_vec = Eigen::VectorXd::Zero(N+1);

  //====================
  // Your code goes here
  //====================

  return rhs_vec;
}  // computeRHS
/* SAM_LISTING_END_3 */

// SOLVE THE LINEAR SYSTEM OF PROBLEM (A)
/* SAM_LISTING_BEGIN_A */
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveA(const Eigen::VectorXd &mesh, FUNCTOR1 &&gamma,
                       FUNCTOR2 &&f) {
  // Nodes are indexed as 0=x_0 < x_1 < ... < x_N = 1
  unsigned N = mesh.size() - 1;
  // Initializations (notice initialization with zeros here)
  Eigen::VectorXd u = Eigen::VectorXd::Zero(N+1);  // solution vec //should be N
  Eigen::SparseMatrix<double> A(N+1, N+1);       // laplacian galerkin mat
  Eigen::SparseMatrix<double> M(N+1, N+1);       // mass galerkin mat
  Eigen::SparseMatrix<double> L(N+1, N+1);       // full galerkin mat

  // I. Build the (full) Galerkin matrix L for the lin. sys.
  //====================
  // Your code goes here
  //====================

//  float h = 1/float(N);
  float h = mesh[1]-mesh[0];
  std::cout<<h<<std::endl;
  std::vector<Eigen::Triplet<double>> A_trip;
  for (int i=0 ; i<N+1 ; ++i){
      A_trip.emplace_back(Eigen::Triplet<double>(i,i,2/h));
  }
  for (int i=0 ; i<N ; ++i){
      A_trip.emplace_back(Eigen::Triplet<double>(i,i+1,-1/h));
  }
  for (int i=1 ; i<N+1 ; ++i){
      A_trip.emplace_back(Eigen::Triplet<double>(i,i-1,-1/h));
  }

  A.setFromTriplets(A_trip.begin(),A_trip.end());
  //std::cout<<A<<std::endl;

  std::vector<Eigen::Triplet<double>> M_trip;
  for (int i=0 ; i<N+1 ; ++i){
    M_trip.emplace_back(Eigen::Triplet<double>(i,i,h*gamma(mesh[i])));
  }
  M.setFromTriplets(M_trip.begin(),M_trip.end());

  L = A+M;
  //std::cout<<L<<std::endl;


  // II. Build the right hand side source vector
  //====================
  // Your code goes here
  //====================

  Eigen::VectorXd v_rhs = Eigen::VectorXd::Zero(N+1);
  for (int i=0 ; i<N+1 ; ++i){
    v_rhs[i] = h*f(mesh[i]);
  }
  int n_tmp = v_rhs.size();
  std::cout<<n_tmp<<std::endl;

  // III. Enforce (zero) dirichlet boundary conditions
  //====================
  // Your code goes here
  //====================

  Eigen::SparseMatrix<double> L_reduced = L.block(1,1,N-1,N-1);
  Eigen::VectorXd v_rhs_reduced = v_rhs.segment(1,N-1);

  std::cout<<"check shapes before LDLT"<<std::endl;
  std::cout<<L_reduced.rows()<<std::endl;
  std::cout<<L_reduced.cols()<<std::endl;
  std::cout<<L_reduced.innerSize()<<std::endl;
  std::cout<<L_reduced.outerSize()<<std::endl;

  std::cout<<"check shapes before LDLT: v_rhs"<<std::endl;
  std::cout<<v_rhs_reduced.size()<<std::endl;


  // IV. Solve the LSE L*u = rhs_vec using an Eigen solver
  //====================
  // Your code goes here
  //====================

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
  solver.compute(L_reduced);
  if (solver.info() != Eigen::Success){
    throw std::runtime_error("Could not decompose the matrix");
  }
  std::cout<<"168"<<std::endl;
  u.segment(1,N-1) = solver.solve(v_rhs_reduced);
  std::cout<<"170"<<std::endl;
  auto ldlt_sol = solver.solve(v_rhs_reduced);
  std::cout<<"check shapes before LDLT: ldlt_sol"<<std::endl;
  std::cout<<ldlt_sol.size()<<std::endl;
  //u.segment(1,N-1) = ldlt_sol;
  u[0] = 666;

  // The solution vector u was initialized with zeros, and therefore already
  // contains the zero Dirichlet boundary data in the first and last entry
  return u;
}  // solveA
/* SAM_LISTING_END_A */

// SOLVE THE LINEAR SYSTEM OF PROBLEM (B)
/* SAM_LISTING_BEGIN_B */
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveB(const Eigen::VectorXd &mesh, FUNCTOR1 &&alpha,
                       FUNCTOR2 &&f, double u0, double u1) {
  // Nodes are indexed as 0=x_0 < x_1 < ... < x_N = 1
  unsigned N = mesh.size() - 1;
  // Initializations
  Eigen::VectorXd u(N+1);                     // solution vec
  Eigen::SparseMatrix<double> A(N+1, N+1);  // laplacian galerkin mat
  // Some tool variables
  double dx_left, dx_right;  // cell widths

  // I. Build the Galerkin matrix A
  //====================
  // Your code goes here
  //====================

  // II. Build the right hand side source vector
  //====================
  // Your code goes here
  //====================

  // III. Enforce dirichlet boundary conditions
  //====================
  // Your code goes here
  //====================

  // IV. Solve the LSE A*u = rhs_vec using an Eigen solver
  //====================
  // Your code goes here
  //====================

  // The solution vector still needs to be supplemented with the known
  // boundary values
  u(0) = u0;  // left boundary node
  u(N) = u1;  // right boundary node
  return u;
}  // solveB
/* SAM_LISTING_END_B */

// Build an sol!ve the LSE corresponding to (C)
/* SAM_LISTING_BEGIN_C */
template <typename FUNCTOR1, typename FUNCTOR2>
Eigen::VectorXd solveC(const Eigen::VectorXd &mesh, FUNCTOR1 &&alpha,
                       FUNCTOR2 &&gamma) {
  // Nodes are indexed as 0=x_0 < x_1 < ... < x_N = 1
  unsigned N = mesh.size() - 1;
  // Initializations (notice initialization with zeros here)
  Eigen::VectorXd u(N+1);                     // solution vec
  Eigen::SparseMatrix<double> A(N+1, N+1);  // laplacian galerkin mat
  Eigen::SparseMatrix<double> M(N+1, N+1);  // mass galerkin mat
  Eigen::SparseMatrix<double> L(N+1, N+1);  // full galerkin mat

  // I. Build the (full) Galerkin matrix L for the lin. sys.
  //====================
  // Your code goes here
  //====================

  // II. Build the right hand side source vector
  //====================
  // Your code goes here
  //====================

  // IV. Solve the LSE A*u = rhs_vec using an Eigen solver
  //====================
  // Your code goes here
  //====================

  return u;
}  // solveC
/* SAM_LISTING_END_C */

}  // namespace LinearFE1D
