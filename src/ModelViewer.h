#pragma once

#include "ModelState.h"
//#include "Gui/DogEditor.h"
#include "DeformationController.h"
#include <igl/opengl/glfw/Viewer.h>

enum ViewMode {
	ViewModeMeshWire = 0,
	ViewModeGauss = 1,
	ViewRulings = 2
};

class ModelViewer {
public:
	ModelViewer(const ModelState& modelState, const DeformationController& DC);

	void render(igl::opengl::glfw::Viewer& viewer);

	ViewMode prevMode;
	ViewMode viewMode;
	bool render_pos_const =  false;

	double rulings_length = 1;
	int rulings_mod = 1;
	double rulings_planar_eps = 0.05;
	bool new_rulings = false;
	bool show_curves = true;
	bool switched_mode;
private:
	void render_mesh_and_wireframe(igl::opengl::glfw::Viewer& viewer);
	void render_gauss_map(igl::opengl::glfw::Viewer& viewer);
	void render_positional_constraints(igl::opengl::glfw::Viewer& viewer);
	void center_and_scale_gauss_sphere(Eigen::MatrixXd& GV, Eigen::MatrixXi& GF);

	void render_mesh(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& Vren, const Eigen::MatrixXi& Fren);
	void reset_state() {first_rendering = true;}
	void clear_edges_and_points(igl::opengl::glfw::Viewer& viewer);

	const ModelState& state;
	const DeformationController& DC;
	bool first_rendering;

	// Used to draw the sphere without directly using opengl routines (but just igl)
	Eigen::MatrixXd sphereV; Eigen::MatrixXi sphereF;
};
