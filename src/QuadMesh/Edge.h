#pragma once
#include <Eigen/Dense>
#include "igl/serialize.h"

#include <vector>

struct Edge  : public igl::Serializable {

  Edge(): v1(-1),v2(-1){}
  Edge(int v1_, int v2_): v1(v1_),v2(v2_) {
    //if (v1 > v2) { std::swap(v1,v2);}
  }
  Edge(const Edge& edge): v1(edge.v1),v2(edge.v2) {}

  void InitSerialization() {
    Add(v1,std::string("v1"));
    Add(v2,std::string("v2"));
  }
  int v1,v2;

  inline bool operator==(const Edge& rhs) const { /* do actual comparison */ 
    if (v1 == rhs.v1 && v2 == rhs.v2) {
      return true;
    }
    if (v1 == rhs.v2 && v2 == rhs.v1) {
      return true;
    }
    return false;
  }

  bool operator<( const Edge& rhs ) const {
    if (v1 < rhs.v1) return true;
    else if (v1 > rhs.v1) return false;
    else return v2 < rhs.v2;
  }
};


struct EdgePoint : public igl::Serializable {
  EdgePoint(){}
  EdgePoint(const Edge& edge, double t) : edge(edge),t(t) {}

  void InitSerialization() {
    Add(edge,std::string("_edge"));
    Add(t,std::string("_edge_t"));
  }

  Eigen::RowVector3d getPositionInMesh(const Eigen::MatrixXd& V) const {return t * V.row(edge.v1) + (1-t) * V.row(edge.v2);}
  Eigen::RowVector3d getPositionInMesh(const Eigen::VectorXd& x, int vn) const {
    Eigen::RowVector3d v1(x(edge.v1),x(vn+edge.v1),x(2*vn+edge.v1));
    Eigen::RowVector3d v2(x(edge.v2),x(vn+edge.v2),x(2*vn+edge.v2));
    return t * v1 + (1-t) * v2;
  }
  static Eigen::MatrixXd getPositionInMesh(const std::vector<EdgePoint>& edgePoints, const Eigen::MatrixXd& V) {
    Eigen::MatrixXd coords(edgePoints.size(),3);
    for (int i = 0; i < edgePoints.size(); i++) {coords.row(i) = edgePoints[i].getPositionInMesh(V);}
      return coords;
  }

  static Eigen::VectorXd getPositionInMesh(const std::vector<EdgePoint>& edgePoints, const Eigen::VectorXd& x) {
    int points_n = edgePoints.size(); Eigen::VectorXd coords(3*points_n);
    for (int i = 0; i < points_n; i++) {
        Eigen::RowVector3d vec = edgePoints[i].getPositionInMesh(x);
        coords(i) = vec(0); coords(points_n+i) = vec(1); coords(2*points_n+i) = vec(2);
    }
    return coords;
  }
  Edge edge;
  double t;
};