#include "Rendering.h"

using namespace std;

void get_wireframe_edges(const Eigen::MatrixXd& V, const QuadTopology& quadTop, Eigen::MatrixXd& E1, Eigen::MatrixXd& E2,
    bool display_border) {

  int e_num = quadTop.E.rows();
  int e_disp_n;
  if (display_border) {
    e_disp_n = e_num;
  } else {
    e_disp_n = 0;
    for (int i = 0; i < e_num; i++) {
      int v1 = quadTop.E(i,0), v2 = quadTop.E(i,1);
      if (!quadTop.is_bnd_v[v1] && !quadTop.is_bnd_v[v2]) {
        e_disp_n++;
      }
    }
  }
  E1.resize(e_disp_n,3); E2.resize(e_disp_n,3);
  int c = 0;
  for (int i = 0; i < e_num; i++) {
    int v1 = quadTop.E(i,0), v2 = quadTop.E(i,1);
    if ((display_border) || (!quadTop.is_bnd_v[v1] && !quadTop.is_bnd_v[v2]) ) {
      E1.row(c) = V.row(v1); E2.row(c) = V.row(v2);
      c++;
    }
  }
}

void render_wireframe_boundary(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const QuadTopology& quadTop,
                        Eigen::RowVector3d color) {
  int e_bnd_num = 0;
  int e_num = quadTop.E.rows();
  for (int i = 0; i < e_num; i++) {
      int v1 = quadTop.E(i,0), v2 = quadTop.E(i,1);
      if (quadTop.is_bnd_v[v1] && quadTop.is_bnd_v[v2]) {
        e_bnd_num++;
      }
  }
  Eigen::MatrixXd E1(e_bnd_num,3), E2(e_bnd_num,3);
  int c = 0;
  for (int i = 0; i < e_num; i++) {
    int v1 = quadTop.E(i,0), v2 = quadTop.E(i,1);
    if (quadTop.is_bnd_v[v1] && quadTop.is_bnd_v[v2]) {
      E1.row(c) = V.row(v1); E2.row(c) = V.row(v2);
      c++;
    }
  }
  viewer.data().add_edges(E1, E2, color);
}

void render_wireframe(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const QuadTopology& quadTop,
    bool display_border) {
  Eigen::MatrixXd E1, E2;
  get_wireframe_edges(V, quadTop, E1, E2, display_border);
  viewer.data().add_edges(E1, E2, Eigen::RowVector3d(0, 0, 0));
}

void plot_vertex_based_rulings(igl::opengl::glfw::Viewer& viewer, const Eigen::MatrixXd& V, const Eigen::MatrixXd& VN, 
    const QuadTopology& quadTop, bool new_rulings, double ruling_length, int rulings_mod, double rulings_planar_eps) {
  Eigen::MatrixXd E1(0,3); Eigen::MatrixXd E2(0,3);
  get_rulings_edges(V, VN, quadTop, new_rulings, ruling_length, rulings_mod, rulings_planar_eps, E1, E2);
  viewer.data().add_edges(E1, E2, Eigen::RowVector3d(0/255,0./255,0./255));
}

void get_rulings_edges(const Eigen::MatrixXd& V, const Eigen::MatrixXd& VN, 
                        const QuadTopology& quadTop, bool new_rulings, double ruling_length, int rulings_mod, double rulings_planar_eps,
                        Eigen::MatrixXd& E1, Eigen::MatrixXd& E2) {
    int r_num = 0;
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
    int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);
    if (!quadTop.is_bnd_v[p_xf_i] && !quadTop.is_bnd_v[p_yf_i] && !quadTop.is_bnd_v[p_xb_i] && !quadTop.is_bnd_v[p_yb_i]) {
      if (p_0_i%rulings_mod == 0) {
        r_num++;
      }
    }
  }

  // plot rulings on non flat areas
  E1.resize(r_num,3);E2.resize(r_num,3);
  int c = 0; //int flat_inner_vertices = 0;
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
    int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);
    if (!quadTop.is_bnd_v[p_xf_i] && !quadTop.is_bnd_v[p_yf_i] && !quadTop.is_bnd_v[p_xb_i] && !quadTop.is_bnd_v[p_yb_i]) {
      if (p_0_i%rulings_mod == 0) {
      
      Eigen::RowVector3d r = get_ruling_direction(VN, p_0_i, p_xf_i, p_xb_i, p_yf_i, p_yb_i, rulings_planar_eps);
      
      double scale = ruling_length*(V.row(p_xf_i)-V.row(p_0_i)).norm();
      E1.row(c) << V(p_0_i,0)-scale*r(0), V(p_0_i,1)-scale*r(1), V(p_0_i,2)-scale*r(2);
      E2.row(c) << V(p_0_i,0)+scale*r(0), V(p_0_i,1)+scale*r(1), V(p_0_i,2)+scale*r(2);
      c++;
      }
    }
  }
}

Eigen::RowVector3d get_ruling_direction(const Eigen::MatrixXd& VN, int p_0_i, int p_xf_i, int p_xb_i, int p_yf_i, int p_yb_i, double rulings_planar_eps) {
   Eigen::RowVector3d n = VN.row(p_0_i);
  Eigen::RowVector3d nx_f = VN.row(p_xf_i);
  Eigen::RowVector3d nx_b = VN.row(p_xb_i);
  Eigen::RowVector3d ny_f = VN.row(p_yf_i);
  Eigen::RowVector3d ny_b = VN.row(p_yb_i);


  // Old rulings
  // make sure the normals are oriented correctly
  if (nx_f.dot(nx_b) < 0) nx_b = -nx_b;
  if (ny_f.dot(ny_b) < 0) ny_b = -ny_b;

  Eigen::RowVector3d nx = nx_f-nx_b;
  Eigen::RowVector3d ny = ny_f-ny_b;

  Eigen::RowVector3d r;

  // if we are flat
  if ( (nx.norm() < rulings_planar_eps) && (ny.norm() < rulings_planar_eps) ) {
    r << 0,0,0; // don't plot any rulings
    //flat_inner_vertices++;
  } else {
    // make sure they are oriented correctly
    if (ny.dot(nx) < 0) ny = -ny;  
    r = n.cross(nx+ny).normalized(); // average cross product and then normalize
  }
  return r;
}
