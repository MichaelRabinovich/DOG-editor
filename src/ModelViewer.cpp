#include "ModelViewer.h"

#include "Gui/Rendering.h"
#include <igl/combine.h>
#include <igl/slice.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>

using namespace std;


ModelViewer::ModelViewer(const ModelState& modelState, const DeformationController& DC) : 
									state(modelState), DC(DC) {
	viewMode = ViewModeMeshWire; 
	prevMode = ViewModeGauss;
	first_rendering = true;

	igl::readOBJ("../data/sphere.obj",sphereV,sphereF);
	center_and_scale_gauss_sphere(sphereV,sphereF);
}

void ModelViewer::render(igl::opengl::glfw::Viewer& viewer) {
	const Dog& dog = DC.getDog();
	clear_edges_and_points(viewer);
	if ((viewMode != prevMode) || (first_rendering)) switched_mode = true;
	prevMode = viewMode;
	if (first_rendering || switched_mode)  {
		viewer.data().clear();
		viewer.core.background_color = Eigen::Vector4f(1, 1, 1, 1);
	}

	if ( viewMode == ViewModeMeshWire ) {
		render_mesh_and_wireframe(viewer);
	} else if (viewMode == ViewModeGauss) {
		render_gauss_map(viewer);
	} else if (viewMode == ViewRulings) {
		render_mesh_and_wireframe(viewer);
		Eigen::MatrixXd VN; //igl::per_vertex_normals(Vren,Fren,VN);
    	igl::per_vertex_normals(dog.getV(),dog.getFTriangular(),VN);
		plot_vertex_based_rulings(viewer, dog.getV(), VN,
						dog.getQuadTopology(), new_rulings, rulings_length, rulings_mod, rulings_planar_eps);
	}

	first_rendering = false;
}

void ModelViewer::clear_edges_and_points(igl::opengl::glfw::Viewer& viewer) {
	viewer.data().set_edges(Eigen::MatrixXd::Zero(0,3), Eigen::MatrixXi::Zero(0,3), Eigen::MatrixXd::Zero(0,3));
	viewer.data().set_points(Eigen::MatrixXd::Zero(0,3), Eigen::MatrixXd::Zero(0,3));
}

void ModelViewer::render_mesh_and_wireframe(igl::opengl::glfw::Viewer& viewer) {
	const Dog* dog = &DC.getDog();
	if (switched_mode) viewer.core.align_camera_center(dog->getVrendering(), dog->getFrendering());
	if (viewMode == ViewModeMeshWire) render_wireframe(viewer, dog->getV(), dog->getQuadTopology());
	render_positional_constraints(viewer);
	render_mesh(viewer, dog->getVrendering(),dog->getFrendering());
}

void ModelViewer::render_mesh(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& Vren, const Eigen::MatrixXi& Fren) {
	if (first_rendering || switched_mode) {
		viewer.data().set_mesh(Vren, Fren);

		Eigen::Vector3d diffuse; diffuse << 0,0,0;
	    Eigen::Vector3d ambient; ambient<< 210.0/255,237.0/255,1.0;
	    Eigen::Vector3d specular; specular << 0,0,0;
	    viewer.data().uniform_colors(ambient,diffuse,specular);
	}
	else {
		 viewer.data().set_vertices(Vren);
    	 viewer.data().compute_normals();
	}

    Eigen::MatrixXd VN; igl::per_vertex_normals(Vren,Fren,VN);
  	viewer.data().set_normals(VN);
}

void ModelViewer::render_positional_constraints(igl::opengl::glfw::Viewer& viewer) {
	const Dog* dog = &DC.getDog();
	Eigen::VectorXd x(dog->getV_vector());
	Eigen::VectorXi b; Eigen::VectorXd bc; DC.get_positional_constraints(b,bc);
	Eigen::VectorXd constrained_pts_coords_vec; igl::slice(x,b,1, constrained_pts_coords_vec);

	int pts_num = b.size()/3;
	Eigen::MatrixXd E1(pts_num,3),E2(pts_num,3);
	for (int i = 0; i < pts_num; i++) {
		E1.row(i) << constrained_pts_coords_vec(i),constrained_pts_coords_vec(pts_num+i),constrained_pts_coords_vec(2*pts_num+i);
		E2.row(i) << bc(i),bc(pts_num+i),bc(2*pts_num+i);
	}
	
	if (render_pos_const) viewer.data().add_edges(E1,E2,Eigen::RowVector3d(1.,0,0));
	//viewer.data().add_points(E1,Eigen::RowVector3d(1.,0,0));
	DC.dogEditor->render_positional_constraints();
}

void ModelViewer::render_gauss_map(igl::opengl::glfw::Viewer& viewer) {
  const Dog* dog = &DC.getDog();
  viewer.data().set_mesh(sphereV, sphereF);
  Eigen::Vector3d diffuse; diffuse << 0.98,0.98,0.98;
  Eigen::Vector3d ambient; ambient << 0,0,0;//0.05*diffuse;
  Eigen::Vector3d specular; specular << 0.05*diffuse;
  viewer.data().uniform_colors(ambient,diffuse,specular);
  //viewer.core.shininess = 0;
  
  if (switched_mode) viewer.core.align_camera_center(dog->getVrendering(), dog->getFrendering());
  //viewer.core.align_camera_center(sphereV, sphereF);
  //viewer.core.show_lines = false;

  Eigen::MatrixXd VN; igl::per_vertex_normals(dog->getVrendering(),dog->getFrendering(),VN);
  //viewer.data.set_normals(VN);
  render_wireframe(viewer,VN,dog->getQuadTopology(), false);
  if (switched_mode) viewer.core.align_camera_center(sphereV, sphereF);
}

void ModelViewer::center_and_scale_gauss_sphere(Eigen::MatrixXd& GV, Eigen::MatrixXi& GF) {
  Eigen::RowVectorXd colmean = GV.colwise().mean();
  for (int i = 0; i < GV.rows(); i++) {
    GV.row(i) = GV.row(i)-colmean;
  }
  Eigen::VectorXd area_v;
  igl::doublearea(GV,GF,area_v);
  double area = area_v.sum()/2.;
  double eps = 2e-1;
  double scale = sqrt((4-eps)*M_PI/area); // make it a little bit smaller so we could see the lines
  GV = GV * scale;
}