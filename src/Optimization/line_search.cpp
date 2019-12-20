#include "line_search.h"

using namespace std;

double line_search(Eigen::VectorXd& x, const Eigen::VectorXd& d, double& step_size, Objective& f, int max_iter) {
	double old_energy = f.obj(x);
  double new_energy = old_energy;
  int cur_iter = 0; int MAX_STEP_SIZE_ITER = min(max_iter,20);
  while (new_energy >= old_energy && cur_iter < MAX_STEP_SIZE_ITER)
  {
    Eigen::VectorXd new_x = x + step_size * d;

    double cur_e = f.obj(new_x);
    
    //cout << "cur_e = " << cur_e << endl;
    if (cur_e >= old_energy)
    {
      step_size /= 2;
      //cout << "step_size = " << step_size << endl;
    }
    else
    {
      x = new_x;
      new_energy = cur_e;
    }
    cur_iter++;
  }
  if (cur_iter <= MAX_STEP_SIZE_ITER) {
  	//cout << "ls success!" << endl;
  } else {
  	cout << "ls failure, is it a local minimum? step size was " << step_size << endl;
  	return old_energy;
  }
  //cout << "step = " << step_size << endl;
  return new_energy;
}

double line_search_l1_directional_derivative(Eigen::VectorXd& x, const Eigen::VectorXd& d, double& step_size, Objective& f, 
      const Constraints& constraints,
        const double& merit_penalty,
        int max_iter) {
  const double c1 = 1e-4;
  ExactL1MeritObjective meritObj(f,constraints,merit_penalty);
  double old_energy = meritObj.obj(x);
  double new_energy = old_energy;
  int cur_iter = 0; int MAX_STEP_SIZE_ITER = min(20,max_iter);
  auto directional_derivative = meritObj.directional_derivative(x,d);

  while ( (new_energy >= old_energy+c1*step_size*directional_derivative) && (cur_iter < MAX_STEP_SIZE_ITER))
  {
    Eigen::VectorXd new_x = x + step_size * d;

    double cur_e = meritObj.obj(new_x);
    
    //cout << "cur_e = " << cur_e << endl;
    if (cur_e >= old_energy)
    {
      step_size *= 0.9;
      //cout << "step_size = " << step_size << endl;
    }
    else
    {
      x = new_x;
      new_energy = cur_e;
    }
    cur_iter++;
  }
  if (cur_iter < MAX_STEP_SIZE_ITER) {
    //cout << "ls success!" << endl;
  } else {
    cout << "ls failure, is it a local minimum?" << endl;
    return old_energy;
  }
  //cout << "step = " << step_size << endl;
  return new_energy;
}



double exact_l2_merit_linesearch(Eigen::VectorXd& x, const Eigen::VectorXd& d, double& step_size, Objective& f, const Constraints& constraints,
				const double& merit_penalty, int max_iter) {

	ExactL2MeritObjective meritObj(f,constraints,merit_penalty);
	return line_search(x,d,step_size, meritObj, max_iter);
}

double exact_l1_merit_linesearch(Eigen::VectorXd& x, const Eigen::VectorXd& d, double& step_size, Objective& f, const Constraints& constraints,
        const double& merit_penalty) {
  ExactL1MeritObjective meritObj(f,constraints,merit_penalty);
  return line_search(x,d,step_size, meritObj); 
}