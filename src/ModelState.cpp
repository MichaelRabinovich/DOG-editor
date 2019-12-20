#include "ModelState.h"

#include <igl/edges.h>
#include <igl/slice.h>
#include <igl/readOBJ.h>
using namespace std;

void ModelState::init_from_mesh(const std::string& mesh_path) {
	std::cout << "Reading mesh " << mesh_path << endl;
	Eigen::MatrixXd V; Eigen::MatrixXi F,F_ren;
	igl::readOBJ(mesh_path, V, F_ren);
	
	// We either read a triangle mesh or a quad mesh
	if (F_ren.cols() == 3 ) {
		F = F_to_Fsqr(F_ren);
	} else {
		F = F_ren;F_ren = Fsqr_to_F(F);
	}
	setup_dog(V,F);
}

void ModelState::setup_dog(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
	QuadTopology quadTop; quad_topology(V,F,quadTop);

	// Scale the mesh so that the entire x curve will be of length 20
	const double edge_l = (V.row(quadTop.bnd_loop[1]) - V.row(quadTop.bnd_loop[0])).norm();
	double cur_v_len = edge_l*sqrt(V.rows());
	auto scaled_V =V* (20. / cur_v_len);

	dog = Dog(scaled_V,F);
}

void ModelState::init_from_planar(int square_h, int square_w) {
	std::cout << "init from planar with  square_h = " << square_h << " square_w = " << square_w << std::endl;
	Eigen::MatrixXd V; Eigen::MatrixXi F;
	get_planar_square_mesh(V, F, square_h, square_w); F = F_to_Fsqr(F);
	setup_dog(V,F);
}