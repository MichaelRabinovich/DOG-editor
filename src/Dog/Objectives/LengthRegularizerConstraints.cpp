#include "LengthRegularizerConstraints.h"

using namespace std;

LengthRegularizerConstraints::LengthRegularizerConstraints(const QuadTopology& quadTop):
		 quadTop(quadTop), vnum(quadTop.v_n) {
	const_n = 2*(quadTop.stars.rows()/5)+1*quadTop.bnd3.rows()/4;
	// Every constraint involves 3 vertices with 9 coordinates
	IJV.resize(12*const_n); // we add the derivatives of the center vertex twice
}

Eigen::VectorXd LengthRegularizerConstraints::Vals(const Eigen::VectorXd& x) const {
	// Edges should be exactly equal
	Eigen::VectorXd constVals(const_n); constVals.setZero();
	
	// Add curve fold constraints
	int const_cnt = 0;
  	#pragma clang loop vectorize(enable)
	//for (auto pair : v_equality_matchings) {
	for (int si = 0; si < quadTop.stars.rows(); si+=5) {
		// get once the x neighbours and once the y neighbours
		for (int coords_i = 0; coords_i < 2; coords_i++) {
				//int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);
			int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1+coords_i), p_xb_i = quadTop.stars(si+3+coords_i);

			const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
			const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
			const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));

			double l1_squared = pow(p0_x-pxf_x,2)+pow(p0_y-pxf_y,2)+pow(p0_z-pxf_z,2);
			double l2_squared = pow(p0_x-pxb_x,2)+pow(p0_y-pxb_y,2)+pow(p0_z-pxb_z,2);
			constVals(const_cnt++) = l1_squared-l2_squared;
		}
	}
	for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
	    int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_xb_i = quadTop.bnd3(si+3);
	    const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
	    const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
	    const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));

	    double l1_squared = pow(p0_x-pxf_x,2)+pow(p0_y-pxf_y,2)+pow(p0_z-pxf_z,2);
		double l2_squared = pow(p0_x-pxb_x,2)+pow(p0_y-pxb_y,2)+pow(p0_z-pxb_z,2);
		constVals(const_cnt++) = l1_squared-l2_squared;
	 }
  if (const_cnt != const_n) {
		cout << "error, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
  }
  //std::cout << "StitchingConstraints avg = " << constVals.cwiseAbs().sum()/constVals.rows() << " and max = " << constVals.cwiseAbs().maxCoeff() << std::endl;
  return constVals;
  
}


void LengthRegularizerConstraints::updateJacobianIJV(const Eigen::VectorXd& x) {

	// Add curve fold constraints
	int const_cnt = 0; int ijv_cnt = 0;
  	#pragma clang loop vectorize(enable)
  	//for (auto pair : v_equality_matchings) {
	for (int si = 0; si < quadTop.stars.rows(); si+=5) {
		// get once the x neighbours and once the y neighbours
		for (int coords_i = 0; coords_i < 2; coords_i++) {
				//int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1), p_yf_i = quadTop.stars(si+2), p_xb_i = quadTop.stars(si+3),p_yb_i = quadTop.stars(si+4);
			int p_0_i = quadTop.stars(si), p_xf_i = quadTop.stars(si+1+coords_i), p_xb_i = quadTop.stars(si+3+coords_i);

			const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
			const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
			const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));

			double l1_squared = pow(p0_x-pxf_x,2)+pow(p0_y-pxf_y,2)+pow(p0_z-pxf_z,2);
			double l2_squared = pow(p0_x-pxb_x,2)+pow(p0_y-pxb_y,2)+pow(p0_z-pxb_z,2);

			double x_der1 = 2*(p0_x-pxf_x), y_der1 = 2*(p0_y-pxf_y), z_der1 = 2*(p0_z-pxf_z);
			double x_der2 = 2*(p0_x-pxb_x), y_der2 = 2*(p0_y-pxb_y), z_der2 = 2*(p0_z-pxb_z);
			
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_0_i,x_der1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_0_i+1*vnum,y_der1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_0_i+2*vnum,z_der1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_xf_i,-x_der1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_xf_i+1*vnum,-y_der1);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_xf_i+2*vnum,-z_der1);

			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_0_i,-x_der2);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_0_i+1*vnum,-y_der2);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_0_i+2*vnum,-z_der2);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_xb_i,+x_der2);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_xb_i+1*vnum,+y_der2);
			IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_xb_i+2*vnum,+z_der2);

			const_cnt++;
		}
	}

	for (int si = 0; si < quadTop.bnd3.rows(); si+=4) {
		// get once the x neighbours and once the y neighbours
		int p_0_i = quadTop.bnd3(si), p_xf_i = quadTop.bnd3(si+1), p_xb_i = quadTop.bnd3(si+3);

		const double p0_x(x(p_0_i+0)); const double p0_y(x(p_0_i+1*vnum)); const double p0_z(x(p_0_i+2*vnum));
		const double pxb_x(x(p_xb_i+0)); const double pxb_y(x(p_xb_i+1*vnum)); const double pxb_z(x(p_xb_i+2*vnum));
		const double pxf_x(x(p_xf_i+0)); const double pxf_y(x(p_xf_i+1*vnum)); const double pxf_z(x(p_xf_i+2*vnum));

		double l1_squared = pow(p0_x-pxf_x,2)+pow(p0_y-pxf_y,2)+pow(p0_z-pxf_z,2);
		double l2_squared = pow(p0_x-pxb_x,2)+pow(p0_y-pxb_y,2)+pow(p0_z-pxb_z,2);

		double x_der1 = 2*(p0_x-pxf_x), y_der1 = 2*(p0_y-pxf_y), z_der1 = 2*(p0_z-pxf_z);
		double x_der2 = 2*(p0_x-pxb_x), y_der2 = 2*(p0_y-pxb_y), z_der2 = 2*(p0_z-pxb_z);
		
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_0_i,x_der1);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_0_i+1*vnum,y_der1);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_0_i+2*vnum,z_der1);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_xf_i,-x_der1);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_xf_i+1*vnum,-y_der1);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_xf_i+2*vnum,-z_der1);

		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_0_i,-x_der2);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_0_i+1*vnum,-y_der2);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_0_i+2*vnum,-z_der2);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_xb_i,+x_der2);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_xb_i+1*vnum,+y_der2);
		IJV[ijv_cnt++] = Eigen::Triplet<double>(const_cnt,p_xb_i+2*vnum,+z_der2);

		const_cnt++;
	}


  if (const_cnt != const_n) {
		cout << "error, const_cnt = " << const_cnt << " but const_n = " << const_n << endl;
		exit(1);
	}
}