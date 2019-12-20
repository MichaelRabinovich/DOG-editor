#pragma once

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/Constraints.h"

#include "../Dog.h"

class LengthRegularizerConstraints: public Constraints {
public:
	LengthRegularizerConstraints(const QuadTopology& quadTop);
	virtual LengthRegularizerConstraints* clone() const {return new LengthRegularizerConstraints(*this);}
	virtual Eigen::VectorXd Vals(const Eigen::VectorXd& x) const;
	virtual void updateJacobianIJV(const Eigen::VectorXd& x);

	virtual void updateLambdaHessianIJV(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda) {
		// Linear constraints have zero second derivative. Empty on purpose
	};

private:
	const QuadTopology& quadTop;
	int vnum;
};
