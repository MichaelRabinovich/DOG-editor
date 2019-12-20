#include "VertexPicker.h"

#include <iostream>
#include <fstream>


#include <igl/project.h>
#include <igl/unproject.h>
#include <igl/point_in_poly.h>
#include <igl/facet_components.h>
#include <igl/barycenter.h>

#include <igl/unproject_in_mesh.h>
#include <igl/Hit.h>

using namespace Eigen;
using namespace igl;
using namespace std;

VertexPicker::VertexPicker(igl::opengl::glfw::Viewer& v,
             const Eigen::MatrixXd &V_,
             const Eigen::MatrixXi &F_tri): V(V_), F(F_tri), viewer(v) {
  // empty on purpose
}


int VertexPicker::pickVertex(int mouse_x, int mouse_y) {
  // Cast a ray in the view direction starting from the mouse position
  double x = mouse_x;
  double y = viewer.core.viewport(3) - mouse_y;
  
  Eigen::RowVector3d pt;
  
  Eigen::Matrix4f modelview = viewer.core.view;// * viewer.core.model;
  int vi = -1;
  
  std::vector<igl::Hit> hits;

  igl::unproject_in_mesh(Eigen::Vector2f(x,y), viewer.core.view,// * viewer.core.model,
      viewer.core.proj, viewer.core.viewport, V, F, pt,hits);

  if (hits.size()> 0) {
    int fi = hits[0].id;
    Eigen::RowVector3d bc;
    bc << 1.0-hits[0].u-hits[0].v, hits[0].u, hits[0].v;
    bc.maxCoeff(&vi);
    vi = F(fi,vi);
  }
  return vi;
}