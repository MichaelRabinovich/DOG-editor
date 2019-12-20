#include "Dog.h"

#include <igl/boundary_loop.h>
#include <igl/edges.h>

using namespace std;
 
Dog::Dog(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) : V(V), F(F), V_ren(V), F_ren(Fsqr_to_F(F)) {
	cout << "init quad top" << endl;
	quad_topology(V,F,quadTop);
	cout << "DOG initialized" << endl;
}

void Dog::update_Vren() {
	V_ren = V;
}