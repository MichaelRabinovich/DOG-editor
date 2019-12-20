#pragma once

#include "../../Optimization/Solver.h"
#include "../../Optimization/Solvers/Pardiso/PardisoSolver.h"

// Equality constrained SQP solver (with a merit function)
class EqSQP : public ConstrainedSolver {
  
public:
	EqSQP(const double& infeasability_epsilon, const double& infeasability_filter, const int& max_newton_iters, const double& merit_p);
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, Constraints& constraints, Eigen::VectorXd& x,
						double convergence_threshold = 1e-4);
private:
	double one_iter(const Eigen::VectorXd& x0, Objective& obj, Constraints& constraints, Eigen::VectorXd& x, double current_merit);
	void build_kkt_system(const Eigen::SparseMatrix<double>& hessian, const Eigen::SparseMatrix<double>& Jacobian,
						Eigen::SparseMatrix<double>& KKT);

	void build_kkt_system_from_ijv(const std::vector<Eigen::Triplet<double> >& hessian_IJV, 
								   const std::vector<Eigen::Triplet<double> >& const_lambda_hessian, int var_n, 
								   const std::vector<Eigen::Triplet<double> >& jacobian_IJV, int const_n);

	double kkt_error(const Eigen::VectorXd& x, Objective& obj, Constraints& eq_constraints);

	const double& infeasability_epsilon; 
	const double& infeasability_filter;
	const int& max_newton_iters;
	const double& merit_p;

	std::vector<Eigen::Triplet<double> > kkt_IJV; Eigen::VectorXi cached_ijv_data;
	Eigen::SparseMatrix<double> A;

  	PardisoSolver m_solver;

	bool first_solve = true;
	Eigen::VectorXd lambda;
	Eigen::VectorXi ai,aj;
	Eigen::VectorXd K;
};