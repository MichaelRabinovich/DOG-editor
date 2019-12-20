#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Objective.h"

// comparing to squared length
class IsometryObjective: public Objective {
  
public:
	IsometryObjective(const QuadTopology& quadTop, const Eigen::VectorXd& x0);
	virtual IsometryObjective* clone() const {return new IsometryObjective(*this);}
	
	virtual double obj(const Eigen::VectorXd& x) const;
	virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const;
	virtual void set_ref(const Eigen::VectorXd& x0);

private:
	virtual void updateHessianIJV(const Eigen::VectorXd& x);
	
	const QuadTopology& quadTop;
	Eigen::VectorXd refL;
	int vnum;
};