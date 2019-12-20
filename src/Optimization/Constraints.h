#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/sparse_cached.h>

class Constraints {
public:

	virtual Constraints* clone() const = 0;

	int getConstNum() const {return const_n;}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const = 0;
	virtual void updateJacobianIJV(const Eigen::VectorXd& x) = 0;
	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) = 0;

	int get_IJV_size() const {return IJV.size();}

	double deviation(const Eigen::VectorXd& x) const {return Vals(x).squaredNorm();}

	const std::vector<Eigen::Triplet<double> >& update_and_get_jacobian_ijv(const Eigen::VectorXd& x) {
    	updateJacobianIJV(x); return IJV; 
  	}
  	const std::vector<Eigen::Triplet<double> >& update_and_get_lambda_hessian(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		updateLambdaHessianIJV(x,lambda); return lambda_hessian_IJV;
  	}

	const std::vector<Eigen::Triplet<double>>& JacobianIJV() {return IJV;}

	// The user just needs to implement JacobianIJV, so other methods could efficiently concatenate multiple constraints jacobian
	virtual const Eigen::SparseMatrix<double>& Jacobian(const Eigen::VectorXd& x) {
		updateJacobianIJV(x);
    	if (cachedJacobian.rows() == 0) {
    		cachedJacobian =  Eigen::SparseMatrix<double>(const_n, x.rows());
      		igl::sparse_cached_precompute(IJV, cached_ijv_data, cachedJacobian);
    	} else {
      		igl::sparse_cached(IJV, cached_ijv_data, cachedJacobian);
    	}
    	return cachedJacobian;
	};

protected:
	int const_n;
	std::vector<Eigen::Triplet<double> > IJV;
	std::vector<Eigen::Triplet<double> > lambda_hessian_IJV;
private:	
	Eigen::SparseMatrix<double> cachedJacobian;
	Eigen::VectorXi cached_ijv_data;
};