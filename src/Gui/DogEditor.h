#pragma once

#include <igl/opengl/glfw/Viewer.h>
#include "VertexPicker.h"

#include "../Dog/Dog.h"

class DogEditor {
public:
	enum EditMode { SELECT_POSITIONAL, TRANSLATE,NONE};//, ROTATE, CUT, GLUE1, GLUE2, NONE };

	DogEditor(igl::opengl::glfw::Viewer& viewer, Dog& dog, EditMode& edit_mode,
         bool& has_new_constraints, Eigen::VectorXi& b, Eigen::VectorXd& bc);
	~DogEditor();

	bool callback_mouse_down();
	bool callback_mouse_move(int mouse_x, int mouse_y);
	bool callback_mouse_up();	

	void render_positional_constraints() const;
	void render_selected_pairs() const;
	void clearHandles();
private:
	Dog& dog;
	EditMode& edit_mode;
	bool& has_new_constraints;
	Eigen::VectorXi& b; Eigen::VectorXd& bc;

	int pick_vertex();
	void select_positional_mouse_down();
	void translate_vertex_edit_mouse_down();

	void onNewHandleID();
	void compute_handle_centroids();
	void get_new_handle_locations(Eigen::Vector3f translation);

	igl::opengl::glfw::Viewer& viewer;
	const Eigen::MatrixXd &V_ren; const Eigen::MatrixXd &V; Eigen::MatrixXd oldV;
	Eigen::MatrixXi F; // The triangular faces of the dog (not dogFrendering)
	VertexPicker vertexPicker;

	//Eigen::Vector3f translation;

	//updated positions of handle vertices, #HV x3
	Eigen::MatrixXd handle_vertex_positions;//(0,3);
	//list of all vertices belonging to handles, #HV x1
	Eigen::VectorXi handle_vertices;//(0,1);
	//for saving constrained vertices
	//vertex-to-handle index, #V x1 (-1 if vertex is free)
	Eigen::VectorXi handle_id;//(0,1);
	int current_handle = -1;

	//index of handle being moved
	int moving_handle = -1;

	//centroids of handle regions, #H x1
	Eigen::MatrixXd handle_centroids;//(0,3);
	int down_mouse_x = -1, down_mouse_y = -1;

	bool action_started = false;
};