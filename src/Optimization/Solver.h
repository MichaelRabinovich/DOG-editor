#pragma once

#include <array>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Objective.h"
#include "Constraints.h"

class Solver {
public:
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve(const Eigen::VectorXd& x0, Objective& obj, Eigen::VectorXd& x) = 0;
};

class ConstrainedSolver {
public:
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, Constraints& constraints, Eigen::VectorXd& x,
			double convergence_threshold = 1e-4) = 0;
};

class IneqConstrainedSolver {
public:
	// x0 is the initial guess, x is the result, the return value is the objective value
	virtual double solve_constrained(const Eigen::VectorXd& x0, Objective& obj, Constraints& eq_constraints,
			Constraints& ineq_constraints, Eigen::VectorXd& x) = 0;
};