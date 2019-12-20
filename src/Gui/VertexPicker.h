#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/Hit.h>

#include <set>

#include "../QuadMesh/Quad.h"

class VertexPicker {
public:
  
  VertexPicker(igl::opengl::glfw::Viewer& v,
         const Eigen::MatrixXd &V_,
         const Eigen::MatrixXi &F_tri);

  int pickVertex(int mouse_x, int mouse_y);

  private:
  const Eigen::MatrixXd &V;
  const Eigen::MatrixXi &F;

  igl::opengl::glfw::Viewer &viewer;
};
