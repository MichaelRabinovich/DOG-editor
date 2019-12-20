#include "DogSolver.h"


using namespace std;

DogSolver::DogSolver(Dog& dog, const Eigen::VectorXd& init_x0, 
        DogSolver::Params& p,
        Eigen::VectorXi& b, Eigen::VectorXd& bc) :

          dog(dog), x(dog.getV_vector()), p(p),
          constraints(dog, init_x0, b, bc),
          obj(dog, init_x0, constraints, p),
          sqpSolver(p.infeasability_epsilon,p.infeasability_filter, p.max_newton_iters, p.merit_p)
           {
    cout << "initializing a dog solver" << endl;
    is_constrained = (b.rows() >0);
}

DogSolver::Constraints::Constraints(const Dog& dog, const Eigen::VectorXd& init_x0,
      Eigen::VectorXi& b, Eigen::VectorXd& bc) : 
                    dogConst(dog.getQuadTopology()),
                    posConst(b,bc),
                    lengthRegConst(dog.getQuadTopology()),
                    compConst({&dogConst}) {
    // Empty on purpose
}

DogSolver::Objectives::Objectives(const Dog& dog, const Eigen::VectorXd& init_x0,
        Constraints& constraints,
          const DogSolver::Params& p) : 
        bending(dog.getQuadTopology(), init_x0), isoObj(dog.getQuadTopology(), init_x0), 
        lengthRegObj(constraints.lengthRegConst, init_x0),
        pointsPosSoftConstraints(constraints.posConst, init_x0),
        compObj(
          {&bending, &isoObj, &lengthRegObj,&pointsPosSoftConstraints},
          {p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(), p.reg_edge_weight/dog.getQuadTopology().E.rows(),p.soft_pos_weight})
          {
    // Empty on purpose
}

void DogSolver::single_iteration(double& constraints_deviation, double& objective) {
  cout << "running a single optimization routine" << endl;
  Eigen::VectorXd x0(x);

  obj.compObj.update_weights({p.bending_weight,p.isometry_weight/dog.getQuadTopology().E.rows(), p.reg_edge_weight/dog.getQuadTopology().E.rows(),p.soft_pos_weight});
  sqpSolver.solve_constrained(x0, obj.compObj, constraints.compConst, x, p.convergence_threshold);

  dog.update_V_vector(x.head(3*dog.get_v_num()));
  
  objective = obj.compObj.obj(x);
  constraints_deviation = constraints.compConst.Vals(x).squaredNorm();
}