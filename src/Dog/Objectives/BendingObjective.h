#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Objective.h"

class BendingObjective: public Objective {
  
public:
	BendingObjective(const QuadTopology& quadTop, const Eigen::VectorXd& x0_init);
	virtual BendingObjective* clone() const {return new BendingObjective(*this);}
	
	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;

private:
	virtual void updateHessianIJV(const Eigen::VectorXd& x);
	
	const QuadTopology& quadTop;
	int vnum;
	std::vector<double> init_edge_lengths;
};