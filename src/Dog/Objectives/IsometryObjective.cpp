#include "IsometryObjective.h"

IsometryObjective::IsometryObjective(const QuadTopology& quadTop, const Eigen::VectorXd& x0)  :
		 quadTop(quadTop), vnum(quadTop.v_n) {
	refL.resize(quadTop.E.rows()); refL.setZero();
	set_ref(x0);
	IJV.resize(quadTop.E.rows()*36);
}

void IsometryObjective::set_ref(const Eigen::VectorXd& x) {

  	int h_cnt = 0;
  	for (int ei = 0; ei < quadTop.E.rows(); ei++) {
		int p_0_i = quadTop.E(ei,0), p_xf_i = quadTop.E(ei,1);

		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));

		double t2 = p0_x-pxf_x;
		double t3 = p0_y-pxf_y;
		double t4 = p0_z-pxf_z;
		refL[h_cnt] = t2*t2+t3*t3+t4*t4; // squared length
		h_cnt++;
  }
}

double IsometryObjective::obj(const Eigen::VectorXd& x) const {
	double e = 0;
	double edge_stretch = 0; double max_stretch = 0;
	int h_cnt = 0;
	for (int ei = 0; ei < quadTop.E.rows(); ei++) {
		int p_0_i = quadTop.E(ei,0), p_xf_i = quadTop.E(ei,1);

		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		double l0 = refL[h_cnt];

  		double t2 = p0_x-pxf_x;
  		double t3 = p0_y-pxf_y;
  		double t4 = p0_z-pxf_z;
  		double t5 = -l0+t2*t2+t3*t3+t4*t4;

  		//double cur_stretch = abs(sqrt(l0)-sqrt(t2*t2+t3*t3+t4*t4))/sqrt(l0);;
  		//edge_stretch += cur_stretch;
  		//max_stretch = std::max(max_stretch,cur_stretch);

  		e += t5*t5;
		h_cnt++;
  }
  //std::cout << "isometry objective = " << e << ", averaged stretch " << edge_stretch/quadTop.E.rows() << ", max stretch = " << max_stretch << std::endl;
  return e;
}

Eigen::VectorXd IsometryObjective::grad(const Eigen::VectorXd& x) const {
  Eigen::VectorXd grad;
  grad.resize(x.rows(),1); grad.setZero();
  int v_num = vnum;
  int h_cnt = 0;

  #pragma clang loop vectorize(enable)
	for (int ei = 0; ei < quadTop.E.rows(); ei++) {
		int p_0_i = quadTop.E(ei,0), p_xf_i = quadTop.E(ei,1);

		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		double l0 = refL[h_cnt];

		double t2 = p0_x-pxf_x;
		double t3 = p0_y-pxf_y;
		double t4 = p0_z-pxf_z;
		double t5 = t2*t2;
		double t6 = t3*t3;
		double t7 = t4*t4;
		double t8 = -l0+t5+t6+t7;
		double t9 = p0_x*2.0;
		double t10 = pxf_x*2.0;
		double t11 = t9-t10;
		double t12 = t8*t11*2.0;
		double t13 = p0_y*2.0;
		double t14 = pxf_y*2.0;
		double t15 = t13-t14;
		double t16 = t8*t15*2.0;
		double t17 = p0_z*2.0;
		double t18 = pxf_z*2.0;
		double t19 = t17-t18;
		double t20 = t8*t19*2.0;

		grad(p_0_i) += t12;
        grad(p_0_i+v_num) += t16;
        grad(p_0_i+2*v_num) += t20;
        grad(p_xf_i) += -t12;
        grad(p_xf_i+v_num) += -t16;
        grad(p_xf_i+2*v_num) += -t20;

        h_cnt++;
  }
  return grad;
}

void IsometryObjective::updateHessianIJV(const Eigen::VectorXd& x) {
  int v_num = vnum;
  int h_cnt = 0;

  int ijv_cnt = 0;
  #pragma clang loop vectorize(enable)
	for (int ei = 0; ei < quadTop.E.rows(); ei++) {
		int p_0_i = quadTop.E(ei,0), p_xf_i = quadTop.E(ei,1);

		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));
		double l0 = refL[h_cnt];

		/*
		Convexication: If l0 <= than current squared length than the hessian is PSD otherwise not,
			and we need to change to a majorizer with an l0* at the hessian which is of a bit of a smaller from l0 (by some epsilon)
		*/
		double x_diff = p0_x-pxf_x;
  		double y_diff = p0_y-pxf_y;
  		double z_diff = p0_z-pxf_z;
  		double cur_squared_l = x_diff*x_diff+y_diff*y_diff+z_diff*z_diff;
  		if (cur_squared_l < l0) {
  			// if l0 is 1 and we are now at 0.9, than by setting l0 to 0.9 we get 0 eigen values
  			//	by setting l0 to 0.9-eps we even get positive ones. 
  			//. This means we need to make l0 smaller by the diff
  			const double eps = 1e-10;
  			l0 -= (l0-cur_squared_l+eps);
  		}

		double t2 = pxf_x*4.0;
		double t3 = pxf_y*4.0;
		double t4 = pxf_z*4.0;
		double t5 = p0_x*4.0;
		double t10 = p0_x*2.0;
		double t11 = pxf_x*2.0;
		double t6 = t10-t11;
		double t7 = p0_x-pxf_x;
		double t8 = p0_y-pxf_y;
		double t9 = p0_z-pxf_z;
		double t12 = t6*t6;
		double t13 = t12*2.0;
		double t14 = t7*t7;
		double t15 = t14*4.0;
		double t16 = t8*t8;
		double t17 = t16*4.0;
		double t18 = t9*t9;
		double t19 = t18*4.0;
		double t20 = p0_y*2.0;
		double t21 = pxf_y*2.0;
		double t22 = t20-t21;
		double t23 = t6*t22*2.0;
		double t24 = p0_z*2.0;
		double t25 = pxf_z*2.0;
		double t26 = t24-t25;
		double t27 = t6*t26*2.0;
		double t28 = p0_y*4.0;
		double t29 = l0*4.0;
		double t30 = t22*t22;
		double t31 = t30*2.0;
		double t32 = t22*t26*2.0;
		double t33 = p0_z*4.0;
		double t34 = t26*t26;
		double t35 = t34*2.0;
		double t36 = -t2+t5;
		double t37 = -t13-t15-t17-t19+t29;
		double t38 = -t3+t28;
		double t39 = -t15-t17-t19+t29-t31;
		double t40 = t15+t17+t19-t29+t31;
		double t41 = -t4+t33;
		double t42 = -t15-t17-t19+t29-t35;
		double t43 = t15+t17+t19-t29+t35;

		// order is p0_x, p0_y, p0_z, pxf_x, pxf_y, pxf_z

		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i,p_0_i, l0*-4.0+t13+t15+t17+t19);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i,p_0_i+vnum, t23);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i,p_0_i+2*vnum, t27);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i,p_xf_i, t37);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i,p_xf_i+vnum, -t23);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i,p_xf_i+2*vnum, -t27);

		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_0_i, t23);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_0_i+vnum, t40);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_0_i+2*vnum, t32);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_xf_i, -t23);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_xf_i+vnum, t39);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i+vnum,p_xf_i+2*vnum, -t32);

		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i, t27);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i+vnum, t32);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_0_i+2*vnum, t43);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i, -t27);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i+vnum, -t32);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_0_i+2*vnum,p_xf_i+2*vnum, t42);

		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i,p_0_i, t37);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i,p_0_i+vnum, -t23);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i,p_0_i+2*vnum, -t27);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i,p_xf_i, t13+t15+t17+t19-t29);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i,p_xf_i+vnum, t23);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i,p_xf_i+2*vnum, t27);

		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_0_i, -t23);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_0_i+vnum, t39);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_0_i+2*vnum, -t32);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i, t23);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i+vnum, t40);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i+vnum,p_xf_i+2*vnum, t32);

		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i, -t27);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i+vnum, -t32);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_0_i+2*vnum, t42);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i, t27);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i+vnum, t32);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(p_xf_i+2*vnum,p_xf_i+2*vnum, t43);

        h_cnt++;
  }
}