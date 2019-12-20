#include "EqSQP.h"

#include "../../QuadMesh/Quad.h"
#include "../../Optimization/line_search.h"

#include "igl/Timer.h"

#include "igl/cat.h"
#include "igl/matlab_format.h"

using namespace std;

EqSQP::EqSQP(const double& infeasability_epsilon, const double& infeasability_filter, const int& max_newton_iters, const double& merit_p) :
    infeasability_epsilon(infeasability_epsilon), infeasability_filter(infeasability_filter), 
             max_newton_iters(max_newton_iters) ,merit_p(merit_p),m_solver(ai,aj,K) 
             {
    cout << "EqSQP::Constructor" << endl;
    m_solver.set_type(-2);
}

double EqSQP::solve_constrained(const Eigen::VectorXd& x0, Objective& f, Constraints& constraints, Eigen::VectorXd& x,
            double convergence_threshold) {
    if (lambda.rows()!= constraints.getConstNum()) {
        lambda.resize(constraints.getConstNum());
        lambda.setZero();
    }

    auto prev_x = x0; double current_merit_p = merit_p; x = x0;
    current_merit_p = max(0.05,lambda.cwiseAbs().maxCoeff()*1.1); // (18.32) in Nocedal
    double ret = f.obj(x0); int iter = 0;

    // smaller tolerance for first iteration (such that typically we will run at least one iteration)
    double error = 10000*kkt_error(x, f, constraints);
    // this might be due to the modified jacobian. We should probably use the normal non modified one for the termination condition..
    while ( (error > infeasability_epsilon) && (iter < max_newton_iters) ) {
        double ret = one_iter(x,f,constraints,x,current_merit_p);
        double const_dev = constraints.Vals(x).norm();
        iter++;
        error = kkt_error(x, f, constraints);
    }
    cout << "finished opt after " << iter << "iterations" << endl;
    return ret;
}

void EqSQP::build_kkt_system_from_ijv(const std::vector<Eigen::Triplet<double> >& hessian_IJV, 
                                     const std::vector<Eigen::Triplet<double> >& const_lambda_hessian, int var_n,
                                     const std::vector<Eigen::Triplet<double> >& jacobian_IJV, int const_n) {
    // build an IJV which contains the hessian, diag matrix along the hessian, jacobian, jacobian transpose,
    //  and diagonal matrix along the J to make sure the whole diagonal is nnz (needed for Pardiso solver)
    int ijv_idx = 0;
    if (kkt_IJV.size() == 0) {
        int kkt_size = hessian_IJV.size() + const_lambda_hessian.size() + 2*jacobian_IJV.size() + var_n + const_n;
        kkt_IJV.resize(kkt_size); int ijv_idx = 0;
    }

    // add hessian+const_lambda_hessian+id*eps_id
    double eps = 1e-8;
    for (int i = 0; i < hessian_IJV.size(); i++) {
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(hessian_IJV[i].row(),hessian_IJV[i].col(),-hessian_IJV[i].value());
    }
    for (int i = 0; i < const_lambda_hessian.size(); i++) {
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(const_lambda_hessian[i].row(),const_lambda_hessian[i].col(),const_lambda_hessian[i].value());
    }
    for (int i = 0; i < var_n; i++) { kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(i,i,-eps);}

    // Add both J and J transpose
    for (int i = 0; i < jacobian_IJV.size(); i++) {
        int row = jacobian_IJV[i].row(), col = jacobian_IJV[i].col(); double val = jacobian_IJV[i].value();
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(col,row+var_n, val); // Jt columed offseted at var_n
        kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(row+var_n, col, val); // J offseted in var_n
    }
    // add zeros along var_n+i,var_n_i
    for (int i = var_n; i < var_n+const_n; i++) kkt_IJV[ijv_idx++] = Eigen::Triplet<double>(i,i,0);

    if ( A.rows() == 0) {
      A =  Eigen::SparseMatrix<double>(var_n+const_n,var_n+const_n);
      igl::sparse_cached_precompute(kkt_IJV, cached_ijv_data, A);
    } else {
      igl::sparse_cached(kkt_IJV, cached_ijv_data, A);
    }
}
void EqSQP::build_kkt_system(const Eigen::SparseMatrix<double>& hessian,
                        const Eigen::SparseMatrix<double>& J, Eigen::SparseMatrix<double>& KKT) {
    Eigen::SparseMatrix<double> id (hessian.rows(),hessian.cols()); id.setIdentity();
    auto eps_id = 1e-8*id;

    Eigen::SparseMatrix<double> H = -hessian - eps_id;
    
    Eigen::SparseMatrix<double> Jt = J.transpose();

    Eigen::SparseMatrix<double> H_jt; igl::cat(2,H,Jt, H_jt);
    Eigen::SparseMatrix<double> zeroM(J.rows(),J.rows());
    Eigen::SparseMatrix<double> J_0; igl::cat(2,J,zeroM,J_0);
    igl::cat(1, H_jt, J_0, KKT);

    //A.makeCompressed();
    Eigen::SparseMatrix<double> id_KKT(KKT.rows(),KKT.rows()); id_KKT.setIdentity(); id_KKT = 0*id_KKT;
    KKT = KKT + id_KKT; // Paradiso's linear solver needs zeros explicitly
}

double EqSQP::one_iter(const Eigen::VectorXd& x0, Objective& f, Constraints& constraints, Eigen::VectorXd& x,
        double current_merit) {
    x = x0;
    int vnum = x.rows()/3;
    double new_e;

    if (lambda.rows()!= constraints.getConstNum()) {
        lambda.resize(constraints.getConstNum());
        lambda.setZero();
    }
    igl::Timer timer; auto init_time = timer.getElapsedTime(); auto t = init_time;
    // Get Hessian
    //auto hessian = f.hessian(x); 
    auto hessian_ijv = f.update_and_get_hessian_ijv(x);
    auto lambda_hessian_ijv = constraints.update_and_get_lambda_hessian(x,lambda);

    // Get Jacobian
    auto jacobian = constraints.Jacobian(x);
    auto jacobian_ijv = constraints.update_and_get_jacobian_ijv(x);

    int var_n = x.rows(); int const_n = constraints.getConstNum();
    build_kkt_system_from_ijv(hessian_ijv, lambda_hessian_ijv, var_n,
                        jacobian_ijv, const_n);

    auto kkt_time = timer.getElapsedTime()-t;  


    double old_e = f.obj(x);
    Eigen::VectorXd g(f.grad(x)); Eigen::VectorXd neg_g = -1*g;
    Eigen::VectorXd d(g.rows()); Eigen::VectorXd lambda_d(lambda.rows());
    if (first_solve) {
        m_solver.iparm[10] = 1; // scaling for highly indefinite symmetric matrices

        m_solver.set_system_matrix(A.triangularView<Eigen::Upper>());
        m_solver.set_pattern();
        m_solver.analyze_pattern();
    
        first_solve = false;
    } else {
        m_solver.update_system_matrix(A.triangularView<Eigen::Upper>());
    }
    m_solver.factorize();

    Eigen::VectorXd constraints_deviation = -1*constraints.Vals(x);
    Eigen::VectorXd rhs_upper = g-jacobian.transpose()*lambda;
    Eigen::VectorXd g_const; igl::cat(1, rhs_upper, constraints_deviation, g_const);
    
    Eigen::VectorXd res;
    m_solver.solve(g_const,res);
    
    for (int d_i = 0; d_i < g.rows(); d_i++) {d[d_i] = res[d_i];}
    for (int d_i = g.rows(); d_i < res.rows(); d_i++) {lambda_d[d_i-g.rows()] = res[d_i];}

    double step_size = 1;
    new_e = exact_l2_merit_linesearch(x,d,step_size,f,constraints,current_merit);
    
    // update lagrange multipliers
    lambda = lambda + step_size*lambda_d;
    auto linesearch_time = timer.getElapsedTime()-t;
    t = timer.getElapsedTime();

    double total_time = timer.getElapsedTime()-init_time;
    return new_e;
}

double EqSQP::kkt_error(const Eigen::VectorXd& x, Objective& obj, Constraints& eq_constraints) {
    auto jacobian = eq_constraints.Jacobian(x);
    
    auto g = obj.grad(x);
    double grad_error = (g-jacobian.transpose()*lambda).norm();
    double eq_const_error = eq_constraints.Vals(x).norm();

    return std::max(grad_error, eq_const_error);
}