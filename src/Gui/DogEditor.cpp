#include "DogEditor.h"

#include "EditingUtils.h"
#include <igl/slice.h>

using namespace std;

DogEditor::DogEditor(igl::opengl::glfw::Viewer& viewer, Dog& dog, EditMode& edit_mode,
         bool& has_new_constraints, Eigen::VectorXi& b, Eigen::VectorXd& bc) : 
		dog(dog), viewer(viewer),
		has_new_constraints(has_new_constraints),b(b),bc(bc),
		V_ren(dog.getVrendering()), V(dog.getV()), F(dog.getFrendering()), 
		vertexPicker(viewer,V_ren,F), edit_mode(edit_mode) {
	
	handle_id.setConstant(V.rows(), 1, -1);	
	oldV = V;
}

bool DogEditor::callback_mouse_down() {
	down_mouse_x = viewer.current_mouse_x;
	down_mouse_y = viewer.current_mouse_y;
	
	if (edit_mode == SELECT_POSITIONAL) {
		select_positional_mouse_down();
	} else if (edit_mode == TRANSLATE) {
		translate_vertex_edit_mouse_down();
	}
	return action_started;
}

bool DogEditor::callback_mouse_move(int mouse_x, int mouse_y) {
	if (!action_started)
		return false;
	
	if (edit_mode == TRANSLATE) {
		Eigen::Vector3f translation = computeTranslation(viewer, mouse_x, down_mouse_x, mouse_y, down_mouse_y, handle_centroids.row(moving_handle));
		get_new_handle_locations(translation);
		down_mouse_x = mouse_x;
		down_mouse_y = mouse_y;
		return true;
	}
	return false;
}

bool DogEditor::callback_mouse_up() {
	if (!action_started)
		return false;
	action_started = false;
	
	if (edit_mode == TRANSLATE) {
		Eigen::Vector3f translation; translation.setZero();
		moving_handle = -1;
		oldV = V;
		
		compute_handle_centroids();
		get_new_handle_locations(translation);
		
		return true;
	}

	if (edit_mode == SELECT_POSITIONAL) {
		return true;
	}
	return false;
};

void DogEditor::onNewHandleID() {
	//store handle vertices too
	int numFree = (handle_id.array()==-1).cast<int>().sum();
	int num_handle_vertices = V.rows() - numFree;
	handle_vertices.setZero(num_handle_vertices);
	handle_vertex_positions.setZero(num_handle_vertices,3);
	
	int count = 0;
	for (long vi = 0; vi<V.rows(); ++vi)
		if(handle_id[vi] >=0)
			handle_vertices[count++] = vi;
	compute_handle_centroids();
	// update handle_vertices_fixed_pos
	igl::slice(V, handle_vertices, 1, handle_vertex_positions);
	// update b and bc (vector form of handle_vertices and handle_vertex_positions)
	const int const_v_num = handle_vertices.rows(); const int v_num = V.rows();

	// b contains normal constraints
	b.resize(3*handle_vertices.rows());
	bc.resize(b.rows());
	for (int i = 0; i < const_v_num; i++) {
		b(i) = handle_vertices(i); bc(i) = handle_vertex_positions(i,0);
		b(const_v_num+i) = handle_vertices(i) + v_num; bc(const_v_num+i) = handle_vertex_positions(i,1);
		b(2*const_v_num+i) = handle_vertices(i) + 2*v_num; bc(2*const_v_num+i) = handle_vertex_positions(i,2);
	}
	has_new_constraints = true;
	cout << "added a new handle" << endl;
}

void DogEditor::get_new_handle_locations(Eigen::Vector3f translation) {
	int count = 0;
	for (long vi = 0; vi<V.rows(); ++vi)
		if(handle_id[vi] >=0) {
			//Eigen::RowVector3f goalPosition = V.row(vi).cast<float>();
			Eigen::RowVector3f goalPosition = oldV.row(vi).cast<float>();
			//Eigen::RowVector3f goalPosition = oldV.row(vi).cast<float>();
			if (handle_id[vi] == moving_handle){
				if( edit_mode == TRANSLATE) goalPosition+=translation;
			}
			handle_vertex_positions.row(count++) = goalPosition.cast<double>();
			oldV.row(vi) = goalPosition.cast<double>();;
		}
	const int const_v_num = handle_vertices.rows(); const int v_num = V.rows();
	bc.resize(b.rows());
	for (int i = 0; i < const_v_num; i++) {
		bc(i) = handle_vertex_positions(i,0);
		bc(const_v_num+i) = handle_vertex_positions(i,1);
		bc(2*const_v_num+i) = handle_vertex_positions(i,2);
	}
}

void DogEditor::compute_handle_centroids() {
	//compute centroids of handles
	int num_handles = handle_id.maxCoeff()+1;
	handle_centroids.setZero(num_handles,3);
	
	Eigen::VectorXi num; num.setZero(num_handles,1);
	for (long vi = 0; vi<V.rows(); ++vi)
	{
		int r = handle_id[vi];
		if ( r!= -1)
		{
			handle_centroids.row(r) += V.row(vi);
			num[r]++;
		}
	}
	
	for (long i = 0; i<num_handles; ++i)
		handle_centroids.row(i) = handle_centroids.row(i).array()/num[i];
	
}

void DogEditor::render_positional_constraints() const {
	Eigen::MatrixXd const_v;
	Eigen::MatrixXd handle_colors(handle_vertex_positions.rows(),3);
    for (int i = 0; i < handle_colors.rows(); i++) {
      if (handle_id[handle_vertices[i]] == moving_handle) {
        handle_colors.row(i) = Eigen::RowVector3d(220./255,0./255,102./255);
      } else {
        handle_colors.row(i) = Eigen::RowVector3d(30./255,80./255,255./255);
      }
    }
    igl::slice(V, handle_vertices, 1, const_v);
    viewer.data().add_points(const_v, handle_colors);
}


void DogEditor::clearHandles() {
	handle_id.setConstant(V.rows(),1,-1);
	handle_vertex_positions.setZero(0,3);
	handle_vertices.resize(0); handle_vertices.setZero(0);

	moving_handle = -1;
	current_handle = -1;
}

void DogEditor::select_positional_mouse_down() {
	int vi = pick_vertex();
	cout << "selected vertex " << vi << endl;
	if (vi >=0) {
		int index = handle_id.maxCoeff()+1;
		if (handle_id[vi] == -1) handle_id[vi] = index;
		current_handle = index;
		
		onNewHandleID();
	}
}

int DogEditor::pick_vertex() {
	int vi = vertexPicker.pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
	if (vi>=0) return vi;
	return -1;
}

void DogEditor::translate_vertex_edit_mouse_down() {
	int vi = pick_vertex();
	if(vi>=0 && handle_id[vi]>=0)  {
		moving_handle = handle_id[vi];
		current_handle = moving_handle;
		oldV = V;
		action_started = true;
	}
}