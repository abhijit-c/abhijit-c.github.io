#include <chrono>
#include <vector>
#include <iostream>

#include <stdlib.h>
#include <math.h>

#include <Eigen/Eigen>

typedef Eigen::Triplet<double> T;
typedef std::chrono::time_point<std::chrono::system_clock> time_pt;

// 1D Poisson equation: u''(x) = f(x), u(0) = u(1) = 0 where f(x) = cos(x).
int main(int argc, char **argv)
{
  if (argc != 2)
  {
    std::cerr << "USAGE: ./1DPoisson.out num_gridpoints\n";
    return -1;
  }
  int N = atoi(argv[1]);
  double h = 1.0/(N-1);

  Eigen::VectorXd f = (M_PI*Eigen::ArrayXd::LinSpaced(N, 0, 1)).sin();
  std::vector<T> coef; coef.reserve(N);

  double sten1 = -1/(h*h), sten2 = 2/(h*h);
  coef.push_back(T(0, 0, sten2)); coef.push_back(T(0,1, sten1));
  for (int k = 1; k < N-1; k++)
  {
    coef.push_back(T(k, k-1, sten1)); coef.push_back(T(k, k+1, sten1));
    coef.push_back(T(k, k, sten2));
  }
  coef.push_back(T(N-1, N-1, sten2)); coef.push_back(T(N-1,N-2, sten1));

  Eigen::SparseMatrix<double> A(N,N);
  A.setFromTriplets(coef.begin(), coef.end());

  time_pt t0 = std::chrono::high_resolution_clock::now();
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
  Eigen::VectorXd u = solver.solve(f);
  time_pt tf = std::chrono::high_resolution_clock::now();

  double error = ( u - (f/(M_PI*M_PI)) ).norm();

  std::cout << "1D Poisson equation where u(0) = u(1) = 0 and f(x) = cos(x).\n" 
            << "Solution computed in "
            << (double)(tf-t0).count()/1e9
            << " s with (2-norm) error "
            << error
            << ".\n"
            << u
            << std::endl;
  return 0;
}
