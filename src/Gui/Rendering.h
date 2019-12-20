#pragma once

#include <igl/opengl/glfw/Viewer.h>
#include "../QuadMesh/Quad.h"
#include "../Dog/Dog.h"

void get_wireframe_edges(const Eigen::MatrixXd& V, const QuadTopology& quadTop, Eigen::MatrixXd& E1, Eigen::MatrixXd& E2, bool display_border);
void render_wireframe(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const QuadTopology& quadTop, bool display_border = true);
void render_wireframe_boundary(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const QuadTopology& quadTop,
		Eigen::RowVector3d color);
void plot_vertex_based_rulings(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const Eigen::MatrixXd& VN,
						const QuadTopology& quad, bool new_rulings, double ruling_length, int rulings_mod = 1, double rulings_planar_eps = 0.05);
void get_rulings_edges(const Eigen::MatrixXd& V, const Eigen::MatrixXd& VN, 
                        const QuadTopology& quadTop, bool new_rulings, double ruling_length, int rulings_mod, double rulings_planar_eps,
                        Eigen::MatrixXd& E1, Eigen::MatrixXd& E2);
Eigen::RowVector3d get_ruling_direction(const Eigen::MatrixXd& VN, int p_0_i, int p_xf_i, int p_xb_i, int p_yf_i, int p_yb_i, double rulings_planar_eps);