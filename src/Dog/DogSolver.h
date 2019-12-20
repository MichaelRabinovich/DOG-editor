#pragma once

#include "Dog.h"

#include "../Optimization/CompositeConstraints.h"
#include "../Optimization/CompositeObjective.h"
#include "../Optimization/PositionalConstraints.h"
#include "../Optimization/QuadraticConstraintsSumObjective.h"
#include "../Optimization/Solvers/EqSQP.h"

#include "Objectives/DogConstraints.h"
#include "Objectives/IsometryObjective.h"
#include "Objectives/LengthRegularizerConstraints.h"
#include "Objectives/BendingObjective.h"

using std::vector;


class DogSolver {
public:
	
	struct Params {
		double bending_weight = 1.;
		double paired_boundary_bending_weight = 1.;
		double isometry_weight = 20000;
		double reg_edge_weight = 0;
		double soft_pos_weight = 5;
		double merit_p = 10;
		int max_newton_iters = 5;
		double infeasability_epsilon = 1e-3;
		double infeasability_filter = 1e-1;
		double convergence_threshold = 1e-6;
	};

	DogSolver(Dog& dog, const Eigen::VectorXd& init_x0, DogSolver::Params& p,
		Eigen::VectorXi& b, Eigen::VectorXd& bc);

	void set_opt_vars(const Eigen::VectorXd& x_i) { x = x_i;}
	Eigen::VectorXd get_opt_vars() { return x;}
	
	void single_iteration(double& constraints_deviation, double& objective);
	void update_point_coords(Eigen::VectorXd& bc) {constraints.posConst.update_coords(bc);}

	Dog& getDog(){return dog;}

	struct Constraints {
		Constraints(const Dog& dog, const Eigen::VectorXd& init_x0, Eigen::VectorXi& b, Eigen::VectorXd& bc);

		DogConstraints dogConst;
		PositionalConstraints posConst;
		LengthRegularizerConstraints lengthRegConst;
		CompositeConstraints compConst;
	};

	struct Objectives {
	  Objectives(const Dog& dog, const Eigen::VectorXd& init_x0,
	  			Constraints& constraints,
	  			const DogSolver::Params& p);

	  	BendingObjective bending;
	  	IsometryObjective isoObj;
	  	QuadraticConstraintsSumObjective lengthRegObj;
      	QuadraticConstraintsSumObjective pointsPosSoftConstraints;
      	CompositeObjective compObj;
	};

private:

	Dog& dog;
	Eigen::VectorXd x; // variables
	bool is_constrained;
	
	// The constraints needs to be defined before the objectives, as some of the objective are dependent on constraints
	DogSolver::Constraints constraints;
	DogSolver::Objectives obj;
	DogSolver::Params& p;

	// Solvers
	EqSQP sqpSolver;
};
