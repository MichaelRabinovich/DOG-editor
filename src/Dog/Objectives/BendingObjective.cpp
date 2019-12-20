#include "BendingObjective.h"

BendingObjective::BendingObjective(const QuadTopology& quadTop, const Eigen::VectorXd& x)  : 
		quadTop(quadTop), vnum(quadTop.v_n) {
	// Number of hessian triplets
	IJV.resize(75*quadTop.stars.rows()/5+27*quadTop.bnd3.rows()/4);
	init_edge_lengths.resize(4*quadTop.stars.rows()/5 + 2*quadTop.bnd3.rows()/4); 
	int cnt = 0;
	for (int si = 0; si < quadTop.stars.rows(); si+=5) {
		int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);

		const double pyb_x(x(p_yb_i+0)); const double pyb_y(x(p_yb_i+1*vnum)); const double pyb_z(x(p_yb_i+2*vnum));
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

		double ex_f_l = sqrt( pow(p0_x-pxf_x,2) + pow(p0_y-pxf_y,2) + pow(p0_z-pxf_z,2) );
		double ex_b_l = sqrt( pow(p0_x-pxb_x,2) + pow(p0_y-pxb_y,2) + pow(p0_z-pxb_z,2) );
		double ey_f_l = sqrt( pow(p0_x-pyf_x,2) + pow(p0_y-pyf_y,2) + pow(p0_z-pyf_z,2) );
		double ey_b_l = sqrt( pow(p0_x-pyb_x,2) + pow(p0_y-pyb_y,2) + pow(p0_z-pyb_z,2) );
		
		// Normalize weights to 1 and choose the ratio as the opposite of the ratio (to have linear precision)
		init_edge_lengths[cnt++] = ex_f_l;
		init_edge_lengths[cnt++] = ex_b_l;
		init_edge_lengths[cnt++] = ey_f_l;
		init_edge_lengths[cnt++] = ey_b_l;
	}

	 for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
	    int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
	    const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
	    const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
	    const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
	    const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

	    double ex_f_l = sqrt( pow(p0_x-pxf_x,2) + pow(p0_y-pxf_y,2) + pow(p0_z-pxf_z,2) );
		double ex_b_l = sqrt( pow(p0_x-pxb_x,2) + pow(p0_y-pxb_y,2) + pow(p0_z-pxb_z,2) );

		init_edge_lengths[cnt++] = ex_f_l;
		init_edge_lengths[cnt++] = ex_b_l;
	}
}

double BendingObjective::obj(const Eigen::VectorXd& x) const {
  double e = 0;
  
  int cnt = 0;  
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
		int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);

		const double pyb_x(x(p_yb_i+0)); const double pyb_y(x(p_yb_i+1*vnum)); const double pyb_z(x(p_yb_i+2*vnum));
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

		double len_ex_f = init_edge_lengths[cnt++];
		double len_ex_b = init_edge_lengths[cnt++];
		double len_ey_f = init_edge_lengths[cnt++];
		double len_ey_b = init_edge_lengths[cnt++];

		double t3 = 1.0/len_ex_b;
		double t4 = 1.0/len_ex_f;
		double t5 = 1.0/len_ey_b;
		double t6 = 1.0/len_ey_f;
		double t2 = t3*(p0_x-pxb_x)+t4*(p0_x-pxf_x)+t5*(p0_x-pyb_x)+t6*(p0_x-pyf_x);
		double t7 = t3*(p0_y-pxb_y)+t4*(p0_y-pxf_y)+t5*(p0_y-pyb_y)+t6*(p0_y-pyf_y);
		double t8 = t3*(p0_z-pxb_z)+t4*(p0_z-pxf_z)+t5*(p0_z-pyb_z)+t6*(p0_z-pyf_z);
		e += t2*t2+t7*t7+t8*t8;


  }

  for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
    int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
    const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
    const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
    const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
    const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

    double len_ex_f = init_edge_lengths[cnt++];
	double len_ex_b = init_edge_lengths[cnt++];
    
	double t3 = 1.0/len_ex_b;
	double t4 = 1.0/len_ex_f;
	double t2 = t3*(p0_x-pxb_x)*2.0+t4*(p0_x-pxf_x)*2.0;
	double t5 = t3*(p0_y-pxb_y)*2.0+t4*(p0_y-pxf_y)*2.0;
	double t6 = t3*(p0_z-pxb_z)*2.0+t4*(p0_z-pxf_z)*2.0;
	e += t2*t2+t5*t5+t6*t6;

  }
  //std::cout << "e = " << std::endl;
  return e;
}

Eigen::VectorXd BendingObjective::grad(const Eigen::VectorXd& x) const {
  Eigen::VectorXd grad;
  grad.resize(x.rows(),1); grad.setZero();
  int v_num = vnum;

  int cnt = 0;
  #pragma clang loop vectorize(enable)
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
        //local_grad.setZero();

        //int p_0_i = i*s+j, p_xf_i = i*s+j+1,p_xb_i = i*s+j-1,p_yf_i = (i+1)*s+j, p_yb_i = (i-1)*s+j;
        int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);
        const double pyb_x(x(p_yb_i+0)); const double pyb_y(x(p_yb_i+1*vnum)); const double pyb_z(x(p_yb_i+2*vnum));
        const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
        const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum));  const double p0_z(x(p_0_i+2*vnum));
        const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
        const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

		double len_ex_f = init_edge_lengths[cnt++];
		double len_ex_b = init_edge_lengths[cnt++];
		double len_ey_f = init_edge_lengths[cnt++];
		double len_ey_b = init_edge_lengths[cnt++];

		double t2 = 1.0/len_ex_b;
		double t3 = 1.0/len_ey_b;
		double t4 = 1.0/len_ex_f;
		double t5 = 1.0/len_ey_f;
		double t6 = t2+t3+t4+t5;
		double t7 = p0_x-pxb_x;
		double t8 = t2*t7;
		double t9 = p0_x-pxf_x;
		double t10 = t4*t9;
		double t11 = p0_x-pyb_x;
		double t12 = t3*t11;
		double t13 = p0_x-pyf_x;
		double t14 = t5*t13;
		double t15 = t8+t10+t12+t14;
		double t16 = p0_y-pxb_y;
		double t17 = t2*t16;
		double t18 = p0_y-pxf_y;
		double t19 = t4*t18;
		double t20 = p0_y-pyb_y;
		double t21 = t3*t20;
		double t22 = p0_y-pyf_y;
		double t23 = t5*t22;
		double t24 = t17+t19+t21+t23;
		double t25 = p0_z-pxb_z;
		double t26 = t2*t25;
		double t27 = p0_z-pxf_z;
		double t28 = t4*t27;
		double t29 = p0_z-pyb_z;
		double t30 = t3*t29;
		double t31 = p0_z-pyf_z;
		double t32 = t5*t31;
		double t33 = t26+t28+t30+t32;

		grad(p_0_i) += t6*t15*2.0;
		grad(p_0_i+v_num) += t6*t24*2.0;
		grad(p_0_i+2*v_num) += t6*t33*2.0;
		grad(p_xb_i+0) += t2*t15*-2.0;
		grad(p_xb_i+v_num) += t2*t24*-2.0;
		grad(p_xb_i+2*v_num) += t2*t33*-2.0;
		grad(p_xf_i) += t4*t15*-2.0;
		grad(p_xf_i+v_num) += t4*t24*-2.0;
		grad(p_xf_i+2*v_num) += t4*t33*-2.0;
		grad(p_yb_i) += t3*t15*-2.0;
		grad(p_yb_i+v_num) += t3*t24*-2.0;
		grad(p_yb_i+2*v_num) += t3*t33*-2.0;
		grad(p_yf_i) += t5*t15*-2.0;
		grad(p_yf_i+v_num) += t5*t24*-2.0;
		grad(p_yf_i+2*v_num) += t5*t33*-2.0;
  }


  for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
	int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
	const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
	const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
	const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
	const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

	double len_ex_f = init_edge_lengths[cnt++];
	double len_ex_b = init_edge_lengths[cnt++];

	double t2 = 1.0/len_ex_b;
	double t3 = 1.0/len_ex_f;
	double t4 = t2*2.0;
	double t5 = t3*2.0;
	double t6 = t4+t5;
	double t7 = p0_x-pxb_x;
	double t8 = t2*t7*2.0;
	double t9 = p0_x-pxf_x;
	double t10 = t3*t9*2.0;
	double t11 = t8+t10;
	double t12 = p0_y-pxb_y;
	double t13 = t2*t12*2.0;
	double t14 = p0_y-pxf_y;
	double t15 = t3*t14*2.0;
	double t16 = t13+t15;
	double t17 = p0_z-pxb_z;
	double t18 = t2*t17*2.0;
	double t19 = p0_z-pxf_z;
	double t20 = t3*t19*2.0;
	double t21 = t18+t20;

	grad(p_0_i) += t6*t11*2.0;
	grad(p_0_i+v_num) += t6*t16*2.0;
	grad(p_0_i+2*v_num) += t6*t21*2.0;
	grad(p_xb_i+0) += t2*t11*-4.0;
	grad(p_xb_i+v_num) += t2*t16*-4.0;
	grad(p_xb_i+2*v_num) += t2*t21*-4.0;
	grad(p_xf_i) += t3*t11*-4.0;
	grad(p_xf_i+v_num) += t3*t16*-4.0;
	grad(p_xf_i+2*v_num) += t3*t21*-4.0;

  }
  return grad;
}


void BendingObjective::updateHessianIJV(const Eigen::VectorXd& x) {
  // Number of ijv values
  int v_num = vnum;
  int ijv_idx = 0; int cnt = 0;
  #pragma clang loop vectorize(enable)
  for (int si = 0; si < quadTop.stars.rows(); si+=5) {
        //local_grad.setZero();

        //int p_0_i = i*s+j, p_xf_i = i*s+j+1,p_xb_i = i*s+j-1,p_yf_i = (i+1)*s+j, p_yb_i = (i-1)*s+j;
        const int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);
        const double pyb_x(x(p_yb_i+0)); const double pyb_y(x(p_yb_i+1*vnum)); const double pyb_z(x(p_yb_i+2*vnum));
        const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
        const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum));  const double p0_z(x(p_0_i+2*vnum));
        const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
        const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

		double len_ex_f = init_edge_lengths[cnt++];
		double len_ex_b = init_edge_lengths[cnt++];
		double len_ey_f = init_edge_lengths[cnt++];
		double len_ey_b = init_edge_lengths[cnt++];

		double t3 = 1.0/len_ex_b;
		double t4 = 1.0/len_ey_b;
		double t5 = 1.0/len_ex_f;
		double t6 = 1.0/len_ey_f;
		double t2 = t3+t4+t5+t6;
		double t7 = t2*t2;
		double t8 = t7*2.0;
		double t9 = 1.0/(len_ex_b*len_ex_b);
		double t10 = t9*2.0;
		double t11 = t3*t5*2.0;
		double t12 = t3*t4*2.0;
		double t13 = t3*t6*2.0;
		double t14 = 1.0/(len_ex_f*len_ex_f);
		double t15 = t14*2.0;
		double t16 = t4*t5*2.0;
		double t17 = t5*t6*2.0;
		double t18 = 1.0/(len_ey_b*len_ey_b);
		double t19 = t18*2.0;
		double t20 = t4*t6*2.0;
		double t21 = 1.0/(len_ey_f*len_ey_f);
		double t22 = t21*2.0;

		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_0_i, t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_xb_i, t2*t3*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_xf_i, t2*t5*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_yb_i, t2*t4*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_yf_i, t2*t6*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_0_i+vnum, t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_xb_i+vnum, t2*t3*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_xf_i+vnum, t2*t5*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_yb_i+vnum, t2*t4*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_yf_i+vnum, t2*t6*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i+2*vnum, t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xb_i+2*vnum, t2*t3*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i+2*vnum, t2*t5*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_yb_i+2*vnum, t2*t4*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_yf_i+2*vnum, t2*t6*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_0_i, t2*t3*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_xb_i, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_xf_i, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_yb_i, t12);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_yf_i, t13);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_0_i+vnum, t2*t3*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xb_i+vnum, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xf_i+vnum, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_yb_i+vnum, t12);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_yf_i+vnum, t13);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_0_i+2*vnum, t2*t3*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xb_i+2*vnum, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xf_i+2*vnum, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_yb_i+2*vnum, t12);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_yf_i+2*vnum, t13);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_0_i, t2*t5*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_xb_i, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_xf_i, t15);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_yb_i, t16);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_yf_i, t17);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_0_i+vnum, t2*t5*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xb_i+vnum, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i+vnum, t15);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_yb_i+vnum, t16);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_yf_i+vnum, t17);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i+2*vnum, t2*t5*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xb_i+2*vnum, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i+2*vnum, t15);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_yb_i+2*vnum, t16);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_yf_i+2*vnum, t17);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i,p_0_i, t2*t4*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i,p_xb_i, t12);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i,p_xf_i, t16);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i,p_yb_i, t19);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i,p_yf_i, t20);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i+vnum,p_0_i+vnum, t2*t4*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i+vnum,p_xb_i+vnum, t12);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i+vnum,p_xf_i+vnum, t16);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i+vnum,p_yb_i+vnum, t19);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i+vnum,p_yf_i+vnum, t20);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i+2*vnum,p_0_i+2*vnum, t2*t4*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i+2*vnum,p_xb_i+2*vnum, t12);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i+2*vnum,p_xf_i+2*vnum, t16);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i+2*vnum,p_yb_i+2*vnum, t19);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yb_i+2*vnum,p_yf_i+2*vnum, t20);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i,p_0_i, t2*t6*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i,p_xb_i, t13);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i,p_xf_i, t17);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i,p_yb_i, t20);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i,p_yf_i, t22);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i+vnum,p_0_i+vnum, t2*t6*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i+vnum,p_xb_i+vnum, t13);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i+vnum,p_xf_i+vnum, t17);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i+vnum,p_yb_i+vnum, t20);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i+vnum,p_yf_i+vnum, t22);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i+2*vnum,p_0_i+2*vnum, t2*t6*-2.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i+2*vnum,p_xb_i+2*vnum, t13);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i+2*vnum,p_xf_i+2*vnum, t17);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i+2*vnum,p_yb_i+2*vnum, t20);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_yf_i+2*vnum,p_yf_i+2*vnum, t22);



   }

   for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
		const int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_yf_i = quadTop.bnd3(si+2), p_xb_i = quadTop.bnd3(si+3);
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		const double pyf_x(x(p_yf_i+0)); const double pyf_y(x(p_yf_i+1*vnum)); const double pyf_z(x(p_yf_i+2*vnum));

		double len_ex_f = init_edge_lengths[cnt++];
		double len_ex_b = init_edge_lengths[cnt++];

		double t3 = 1.0/len_ex_b;
		double t4 = t3*2.0;
		double t5 = 1.0/len_ex_f;
		double t6 = t5*2.0;
		double t2 = t4+t6;
		double t7 = t2*t2;
		double t8 = t7*2.0;
		double t9 = 1.0/(len_ex_b*len_ex_b);
		double t10 = t9*8.0;
		double t11 = t3*t5*8.0;
		double t12 = 1.0/(len_ex_f*len_ex_f);
		double t13 = t12*8.0;

		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_0_i, t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_xb_i, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i,p_xf_i, t2*t5*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_0_i+vnum, t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_xb_i+vnum, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+vnum,p_xf_i+vnum, t2*t5*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i+2*vnum, t8);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xb_i+2*vnum, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i+2*vnum, t2*t5*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_0_i, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_xb_i, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i,p_xf_i, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_0_i+vnum, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xb_i+vnum, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+vnum,p_xf_i+vnum, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_0_i+2*vnum, t2*t3*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xb_i+2*vnum, t10);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xb_i+2*vnum,p_xf_i+2*vnum, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_0_i, t2*t5*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_xb_i, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i,p_xf_i, t13);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_0_i+vnum, t2*t5*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xb_i+vnum, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i+vnum, t13);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i+2*vnum, t2*t5*-4.0);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xb_i+2*vnum, t11);
		IJV[ijv_idx++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i+2*vnum, t13);
	}
}