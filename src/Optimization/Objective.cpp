#include "Objective.h"

using namespace std;

// builds the hessian from an IJV
const Eigen::SparseMatrix<double>& Objective::hessian(const Eigen::VectorXd& x) {
  updateHessianIJV(x);
  if (!is_H_cached) {
    cachedH =  Eigen::SparseMatrix<double>(x.rows(),x.rows());
    igl::sparse_cached_precompute(IJV, cached_ijv_data, cachedH);
    is_H_cached = true;
  } else {
    igl::sparse_cached(IJV, cached_ijv_data, cachedH);
  }
  return cachedH;
}

std::vector<Eigen::Triplet<double>> Objective::to_triplets(Eigen::SparseMatrix<double> & M) {
    std::vector<Eigen::Triplet<double>> v;
    for(int i = 0; i < M.outerSize(); i++)
        for(typename Eigen::SparseMatrix<double>::InnerIterator it(M,i); it; ++it)
            v.emplace_back(it.row(),it.col(),it.value());
    return v;
  }

void Objective::check_grad(const Eigen::VectorXd& x) const {
  	Eigen::VectorXd g, fin_g;
  	g = grad(x);
  	
    finiteGradient(x, fin_g);
    std::cout << "Checking grad!" << std::endl;
    std::cout << "g.norm() = " << g.norm() << std::endl;
    std::cout << "fin_g.norm() = " << fin_g.norm() << std::endl;
    std::cout << "(g-fin_g).norm() = " << (g-fin_g).norm() << std::endl;
    //std::cout << "g.rows() = " << g.rows() << std::endl;
    //std::cout << "fin_g.rows() = " << fin_g.rows() << std::endl;
    
	/*
    for (int i = 0; i < g.rows(); i++) {
    	if (abs(g(i)) > 1e-12) {
    		cout << "g(" << i << ") = " << g(i) << endl;
    	}
    	if (abs(fin_g(i)) > 1e-12) {
    		cout << "fin_g(" << i << ") = " << fin_g(i) << endl;
    	}
    }
    //exit(1);
	*/
  }

void Objective::finiteGradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad, int accuracy) const {
    // accuracy can be 0, 1, 2, 3
    const double eps = 2.2204e-6;
    static const std::array<std::vector<double>, 4> coeff =
    { { {1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3} } };
    static const std::array<std::vector<double>, 4> coeff2 =
    { { {1, -1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4} } };
    static const std::array<double, 4> dd = {2, 12, 60, 840};

    grad.resize(x.rows());grad.setZero();
    Eigen::VectorXd& xx = const_cast<Eigen::VectorXd&>(x);

    const int innerSteps = 2*(accuracy+1);
    const double ddVal = dd[accuracy]*eps;

    for (int d = 0; d < x.rows(); d++) {
      grad[d] = 0;
      for (int s = 0; s < innerSteps; ++s)
      {
        double tmp = xx[d];
        xx[d] += coeff2[accuracy][s]*eps;
        grad[d] += coeff[accuracy][s]*obj(xx);
        xx[d] = tmp;
      }
      grad[d] /= ddVal;
    }
  }