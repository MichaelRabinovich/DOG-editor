#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/slice.h>

#include "Constraints.h"

class PositionalConstraints : public Constraints {
public:
	PositionalConstraints(const Eigen::VectorXi& b, const Eigen::VectorXd& bc) : b(b),bc(bc) {
		const_n = b.rows();
		IJV.resize(const_n);
	};

	virtual PositionalConstraints* clone() const {return new PositionalConstraints(*this);}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const {
		Eigen::VectorXd constrained_pts_coords; igl::slice(x,b,1, constrained_pts_coords);
		// A vector of x(b)-bc for the constraints x(b)-bc = 0
		return constrained_pts_coords-bc;
	}
	virtual void updateJacobianIJV(const Eigen::VectorXd& x) {
		int const_n = 0;
		for (int b_i = 0; b_i < b.rows(); b_i++ ) {
			int var_const_idx = b(b_i);
			// Set the derivative at the 'var_const_idx' as d(x(val_idx)-value)/d(val_idx) = 1
			IJV[b_i] = Eigen::Triplet<double>(const_n++, var_const_idx, 1);
		}
	}
	
	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// Linear constraints have zero second derivative. Empty on purpose
	};

	void update_coords(Eigen::VectorXd& bc_i) {bc = bc_i;}
	Eigen::VectorXi getPositionIndices() const {return b;}
	Eigen::VectorXd getPositionVals() const {return bc;}

private:
	Eigen::VectorXi b; Eigen::VectorXd bc;
};
/*
// Not done with QuadraticConstraintsSumObjective over PositionalConstraints because we need the hessian here
// and it was not implemented (with the chain rule) for QuadraticConstraintsSumObjective yet
class SoftPositionalConstraints : public Objective {
public:
	SoftPositionalConstraints(const Eigen::VectorXi& b, const Eigen::VectorXd& bc) : posConst(b,bc), innerObj(posConst) {}

	double obj(const Eigen::VectorXd& x) const { return innerObj.obj(x);}
	
	Eigen::VectorXd grad(const Eigen::VectorXd& x) const {return innerObj.grad(x);}
	virtual Eigen::SparseMatrix<double> hessian(const Eigen::VectorXd& x) const {
		// The constraints is linear so the hessian is just simple 2*J'*J
		// (Here it's even simpler as just 2*(diagonals with one at constrained indices) but this is simpler to write)
		auto J = posConst.Jacobian(x);
		return 2*J.transpose()*J;
	}
private:
	const PositionalConstraints posConst;
};*/