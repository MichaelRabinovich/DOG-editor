#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <Eigen/Dense>
#include <vector>

#include "igl/serialize.h"

#include "../QuadMesh/Quad.h"

// Encapsulate a dog mesh and its rendered mesh
class Dog {
public:
	Dog(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
	Dog(){}

	void update_Vren();
	void update_V(const Eigen::MatrixXd& V_new) {V = V_new; update_Vren();}
	void update_V_vector(const Eigen::VectorXd& x) {vec_to_mat2(x,V); update_Vren();}

	Dog* get_submesh(int submesh_i);

	int get_v_num() {return V.rows();}

	const Eigen::MatrixXi& getF() const {return F;}
	const Eigen::MatrixXd& getV() const {return V;}
	//const Eigen::MatrixXd& getFlatV() const {return flatV;} somehow not working now
	const QuadTopology& getQuadTopology() const {return quadTop;}
	Eigen::MatrixXd& getVMutable() {return V;}
	Eigen::VectorXd getV_vector() const {Eigen::VectorXd x; mat2_to_vec(V,x); return x;}
	Eigen::MatrixXi getFTriangular() const {return Fsqr_to_F(F);} // useful for the editor who needs a triangular mesh (still different then the rendering)
	const Eigen::MatrixXi& getFrendering() const {return F_ren;}
	const Eigen::MatrixXd& getVrendering() const {return V_ren;}

private:
	// The quad mesh
	Eigen::MatrixXd V; Eigen::MatrixXi F;
	// Indices of boundary curves (also when there are creases). Only relevant for square patches and used for wallpapers.
	QuadTopology quadTop;
	// The initial rendered (triangular) mesh
	Eigen::MatrixXd V_ren; Eigen::MatrixXi F_ren;
};
