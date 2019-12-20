#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Constraints.h"

/*
A constraint set built from a set of other constraints
*/
class CompositeConstraints : public Constraints {
public:
	CompositeConstraints() {};
	virtual CompositeConstraints* clone() const {return new CompositeConstraints(*this);}

	CompositeConstraints(const std::vector<Constraints*>& constraints_i, const int var_range = -1) : var_range(var_range){
		constraints.resize(constraints_i.size());
		for (int i = 0; i < constraints.size(); i++) constraints[i] = constraints_i[i];
		const_n = 0; 
		for (auto cnst: constraints) {const_n+=cnst->getConstNum(); ijv_size += cnst->get_IJV_size();}
		IJV.resize(ijv_size);
	};

	void add_constraints_permanent(Constraints* cnst) {
		constraints.push_back(cnst);
		const_n+= cnst->getConstNum();
		ijv_size += cnst->get_IJV_size();
		IJV.resize(ijv_size);
	}

	void add_constraints(Constraints* cnst) {
		constraints.push_back(cnst->clone());
		const_n+= cnst->getConstNum();
		ijv_size += cnst->get_IJV_size();
		IJV.resize(ijv_size);
	}

	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x_whole) const {
		Eigen::VectorXd x;
		if (var_range == -1) x = x_whole; else x = x_whole.head(var_range);
		Eigen::VectorXd vals(const_n);
		int const_cnt = 0; 
		for (auto cnst: constraints) {
			auto cnst_vals = cnst->Vals(x);
			for (int val_const_i = 0; val_const_i < cnst_vals.rows(); val_const_i++) {vals[const_cnt++] = cnst_vals[val_const_i];}
		}
		if (const_n!= const_cnt) {
			std::cout << "Error! const_n = " << const_n << " and const_cnt = " << const_cnt << std::endl;
			exit(1);
		}
		return vals;
	}
	virtual void updateJacobianIJV(const Eigen::VectorXd& x_whole) {
		Eigen::VectorXd x;
		if (var_range == -1) x = x_whole; else x = x_whole.head(var_range);
		int row_base = 0; int ijv_idx = 0;
		for (auto cnst: constraints) {
			cnst->updateJacobianIJV(x);
			const std::vector<Eigen::Triplet<double> >& cnst_IJV = cnst->JacobianIJV();
			
			for (auto val : cnst_IJV) {
				IJV[ijv_idx++] = Eigen::Triplet<double>(val.row() + row_base, val.col(), val.value());
			}
			row_base += cnst->getConstNum();
		}
	}
	
	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x_whole, const Eigen::VectorXd& lambda){
		Eigen::VectorXd x;
		if (var_range == -1) x = x_whole; else x = x_whole.head(var_range);
		int lambda_idx = 0;
		// Each set of constraints has a different sub set of the lagrange multipliers
		for (auto cnst: constraints) {
			// Get the relevant lagrange multipliers for the subset of constraints
			int lambda_size = cnst->getConstNum(); Eigen::VectorXd cnst_lambda(lambda_size);
			for (int i = 0; i < lambda_size; i++) { cnst_lambda(i) = lambda(lambda_idx+i);}
			// Update the lambda hessian ijv of this constraint
			cnst->updateLambdaHessianIJV(x,cnst_lambda);
			// set the next abse index for the lagrange multipliers
			lambda_idx += lambda_size;
		}
	}
private:
	std::vector<Constraints*> constraints;
	int var_range;
	int ijv_size = 0;
};