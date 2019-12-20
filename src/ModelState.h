#pragma once

#include "Dog/Dog.h"
#include "DeformationController.h"

#include "igl/serialize.h"

struct ModelState {
	Dog dog;
	DeformationController DC;

	void init_from_mesh(const std::string& mesh_path);
	void init_from_planar(int square_h, int square_w);

private:
	void setup_dog(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
};