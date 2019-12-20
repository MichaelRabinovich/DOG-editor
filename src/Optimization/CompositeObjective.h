#pragma once

#include "Objective.h"
#include <igl/Timer.h>

/*
An objective built from a weighted composition of other objectives
*/
class CompositeObjective: public Objective {
  
public:
	CompositeObjective(){/*empty on purpose*/}
	CompositeObjective(const std::vector<Objective*>& objectives_i, const std::vector<double>& weights, const int var_range = -1) : 
				weights(weights), var_range(var_range) {
		objectives.resize(objectives_i.size());
		for (int i = 0; i < objectives.size(); i++) {
			objectives[i] = objectives_i[i]->clone();
			ijv_size += objectives[i]->get_hessian_IJV_size();
		}
		IJV.resize(ijv_size);
	}

	void update_weights(const std::vector<double>& weights_i) { 
		weights = weights_i; 
		if (weights.size()<objectives.size()) {
			std::cout << "Error in CompositeObjective::update_weights - " << "There are " << objectives.size() << " objectives" << 
				" but got " << weights.size() << " new weights " << std::endl;

		}}

	virtual CompositeObjective* clone() const {return new CompositeObjective(*this);}
	void add_objective(Objective* e, double w = 1.) {
		objectives.push_back(e->clone());
		weights.push_back(w);
		std::cout<<std::endl;
		ijv_size += e->get_hessian_IJV_size();
		IJV.resize(ijv_size);
	}

	void add_objective_permanent(Objective& e, double w = 1.) {
		objectives.push_back(&e); 
		weights.push_back(w);
		ijv_size += e.get_hessian_IJV_size();
		IJV.resize(ijv_size);
	}

	virtual double obj(const Eigen::VectorXd& x_whole) const {
		Eigen::VectorXd x;
		if (var_range == -1) x = x_whole; else x = x_whole.head(var_range);
		double obj = 0;
		for (int i = 0; i < objectives.size(); i++) {
			if (weights[i]>0) obj+=weights[i]*objectives[i]->obj(x);
		}
		return obj;
	}

	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x_whole) const {
		Eigen::VectorXd x;
		if (var_range == -1) x = x_whole; else x = x_whole.head(var_range);
		Eigen::VectorXd grad(x.rows()); grad.setZero();
		for (int i = 0; i < objectives.size(); i++) {
			if (weights[i]>0) grad+=weights[i]*objectives[i]->grad(x);
		}
		if (var_range == -1) {
			return grad;
		} else {
			Eigen::VectorXd whole_grad(x_whole.rows()); whole_grad.setZero();
			whole_grad.head(var_range) = grad;
			return whole_grad;
		}
	};

	virtual void set_ref(const Eigen::VectorXd& x0_whole) {
		Eigen::VectorXd x0;
		if (var_range == -1) x0 = x0_whole; else x0 = x0_whole.head(var_range);
		for (auto obj : objectives) {obj->set_ref(x0);}
	};


private:
	virtual void updateHessianIJV(const Eigen::VectorXd& x_whole) {
		Eigen::VectorXd x;
		if (var_range == -1) x = x_whole; else x = x_whole.head(var_range);
		int ijv_idx = 0;
		for (int i = 0; i < objectives.size(); i++) {
			const std::vector<Eigen::Triplet<double> >& obj_IJV = objectives[i]->update_and_get_hessian_ijv(x);
			for (auto val : obj_IJV) {IJV[ijv_idx++] = Eigen::Triplet<double>(val.row(),val.col(),weights[i]*val.value());}
		}
	}

	std::vector<Objective*> objectives;
	std::vector<double> weights;
	int var_range;
	int ijv_size = 0;
};