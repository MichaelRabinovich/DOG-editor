#pragma once

#include "igl/serialize.h"

#include "Dog/DogSolver.h"
#include "Gui/DogEditor.h"

class DeformationController {
public:
	DeformationController();
	~DeformationController() {if (dogSolver) delete dogSolver;}
	void init_from_new_dog(Dog& dog);
	void init_viewer(igl::opengl::glfw::Viewer& viewer_i) {viewer = &viewer_i;}	

	bool has_constraints();

	void get_positional_constraints(Eigen::VectorXi& b_out, Eigen::VectorXd& bc_out) const {b_out=b;bc_out = bc;};
	void single_optimization();
	void reset_constraints();
	const Dog& getDog() const {return dogSolver->getDog();}

	DogEditor::EditMode edit_mode = DogEditor::NONE;
	bool has_new_constraints = false;
	
	DogEditor* dogEditor;
	DogSolver::Params p;
	double constraints_deviation;
	double objective;
private:
	void reset_dog_solver();

	void update_point_coords(Eigen::VectorXd& bc_i);
	void add_positional_constraints(const Eigen::VectorXi& new_b, const Eigen::VectorXd& new_bc);

	igl::opengl::glfw::Viewer* viewer;

	// Points positional constraints
	Eigen::VectorXi b; Eigen::VectorXd bc;
	Eigen::VectorXd init_x0;
	// This needs to reset sometimes. 
	// For instance when a new soft constraint is added (but not when the constraint value change), or when a entirely new DOG is loaded
	//	Since this amounts to a different objective/hessian sparsity pattern
	DogSolver* dogSolver;
};