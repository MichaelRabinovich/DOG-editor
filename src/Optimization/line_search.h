#pragma once

#include "Objective.h"
#include "Constraints.h"
//#include <Eigen/Array.h>

double line_search(Eigen::VectorXd& x, const Eigen::VectorXd& d, double& step_size, Objective& f, int max_iter = 20);

double line_search_l1_directional_derivative(Eigen::VectorXd& x, const Eigen::VectorXd& d, double& step_size, Objective& f, 
      const Constraints& constraints,
        const double& merit_penalty,
        int max_iter = 20);

double exact_l2_merit_linesearch(Eigen::VectorXd& x, const Eigen::VectorXd& d, double& step_size, Objective& f, const Constraints& constraints,
				const double& merit_penalty, int max_iter = 20);

double exact_l1_merit_linesearch(Eigen::VectorXd& x, const Eigen::VectorXd& d, double& step_size, Objective& f, const Constraints& constraints,
				const double& merit_penalty);


class ExactL1MeritObjective: public Objective {
  
public:
	ExactL1MeritObjective(const Objective& obj, const Constraints& constraints, const double& merit_p)  : inner_obj(obj), inner_constraints(constraints), merit_p(merit_p) {};
	virtual ExactL1MeritObjective* clone() const {return new ExactL1MeritObjective(*this);}
	
	virtual double obj(const Eigen::VectorXd& x) const {
		return inner_obj.obj(x)+merit_p*inner_constraints.Vals(x).lpNorm<1>();
	};
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const {
		// Should not get here..
		Eigen::VectorXd g(x.rows());
		std::cout << "Error, L2Merit function is not differentiable" << std::endl; exit(1);
		return g;
	}

	double directional_derivative(const Eigen::VectorXd& x, const Eigen::VectorXd& d) {
		return inner_obj.grad(x).dot(d)-merit_p*inner_obj.obj(x);
	}

private:
	const Objective& inner_obj; 
	const Constraints& inner_constraints;
	const double& merit_p;
};

class ExactL2MeritObjective: public Objective {
  
public:
	ExactL2MeritObjective(const Objective& obj, const Constraints& constraints, const double& merit_p)  : inner_obj(obj), inner_constraints(constraints), merit_p(merit_p) {};
	virtual ExactL2MeritObjective* clone() const {return new ExactL2MeritObjective(*this);}
	
	virtual double obj(const Eigen::VectorXd& x) const {
		return inner_obj.obj(x)+merit_p*inner_constraints.Vals(x).norm();
	};
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const {
		// Should not get here..
		Eigen::VectorXd g(x.rows());
		std::cout << "Error, L2Merit function is not differentiable" << std::endl; exit(1);
		return g;
	}

private:
	const Objective& inner_obj; 
	const Constraints& inner_constraints;
	const double& merit_p;
};