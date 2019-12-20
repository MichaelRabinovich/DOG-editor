#pragma once

#include <array>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/sparse_cached.h>


/*
A class for objectives. Base classes need to implement obj(), grad() and clone() (see QuadraticConstraintsSumObjective for an example)
*/
class Objective {
  
public:
  virtual Objective* clone() const = 0;

  virtual double obj(const Eigen::VectorXd& x) const = 0;
  virtual Eigen::VectorXd grad(const Eigen::VectorXd& x) const = 0;

  int get_hessian_IJV_size() {return IJV.size();}
  
  const std::vector<Eigen::Triplet<double> >& update_and_get_hessian_ijv(const Eigen::VectorXd& x) {
    updateHessianIJV(x); return IJV; 
  }
  // For objectives that depeends on a state/reference shape (which might also be updated at various times)
  virtual void set_ref(const Eigen::VectorXd& x0) {};

  void check_grad(const Eigen::VectorXd& x) const;

  std::vector<Eigen::Triplet<double>> to_triplets(Eigen::SparseMatrix<double> & M);
  // builds the hessian from an IJV
  virtual const Eigen::SparseMatrix<double>& hessian(const Eigen::VectorXd& x);
 protected:
   std::vector<Eigen::Triplet<double> > IJV; Eigen::VectorXi cached_ijv_data;
   Eigen::SparseMatrix<double> cachedH;
   bool is_H_cached = false;
 private:

  
  // return 0 sparse matrix if not implemented
  virtual void updateHessianIJV(const Eigen::VectorXd& x) { /*empty on purpose */ } 

  void finiteGradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad, int accuracy = 1) const;
};