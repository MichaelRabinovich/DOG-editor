#include "DeformationController.h"

#include <queue>
using namespace std;

DeformationController::DeformationController() : dogEditor(NULL), dogSolver(NULL) {
	// empty on piurpose
}

void DeformationController::init_from_new_dog(Dog& dog) {
	cout << "init_from_new_dog" << std::endl;
	dogEditor = new DogEditor(*viewer, dog, edit_mode, has_new_constraints,b,bc);

	init_x0 = dog.getV_vector();
	if (dogSolver) delete dogSolver;
	dogSolver = new DogSolver(dog,init_x0, p, b, bc);
}

bool DeformationController::has_constraints() {
	return b.rows();
}

void DeformationController::single_optimization() {
	if (has_new_constraints) reset_dog_solver();
	dogSolver->update_point_coords(bc);
	dogSolver->single_iteration(constraints_deviation, objective);
}

void DeformationController::update_point_coords(Eigen::VectorXd& bc_i) {
	bc = bc_i; dogSolver->update_point_coords(bc);
}

void DeformationController::add_positional_constraints(const Eigen::VectorXi& new_b, const Eigen::VectorXd& new_bc) {
	Eigen::VectorXi old_b = b; Eigen::VectorXd old_bc = bc;
	b.resize(old_b.rows()+new_b.rows()); bc.resize(b.rows());
	if (old_b.rows()) {
		b << old_b,new_b; bc << old_bc,new_bc;
	} else {
		b = new_b; bc = new_bc; // Eigen's concatenate crashes if one of them is empty
	}
	has_new_constraints = true;
}

void DeformationController::reset_constraints() {
	b.resize(0);
	bc.resize(0); 
	dogEditor->clearHandles(); 
	reset_dog_solver();
}

void DeformationController::reset_dog_solver() {
	Dog& dog = dogSolver->getDog();
	auto vars = dogSolver->get_opt_vars();
	if (dogSolver) delete dogSolver;
	cout << "resetting dog solver" << endl;
	dogSolver = new DogSolver(dog,init_x0, p, b, bc);
	dogSolver->set_opt_vars(vars);
	has_new_constraints = false;
}