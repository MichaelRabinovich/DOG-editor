#pragma once

#include <igl/opengl/glfw/Viewer.h>

Eigen::Vector3f computeTranslation(igl::opengl::glfw::Viewer& viewer, int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);
Eigen::Vector4f computeRotation(igl::opengl::glfw::Viewer& viewer, int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);