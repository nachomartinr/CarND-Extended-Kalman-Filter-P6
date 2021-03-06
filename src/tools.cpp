#include <iostream>
#include "tools.h"

using namespace std;

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
	TODO:
		* Calculate the RMSE here.
	*/
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of inputs
  if((estimations.size() != ground_truth.size()) || (estimations.size() == 0)){
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  //accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){

    VectorXd residual = estimations[i] - ground_truth[i];

    //square coefficient-wise
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse/estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
	TODO:
  * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4);
  //recover state parameters
  const float px = x_state(0);
  const float py = x_state(1);
  const float vx = x_state(2);
  const float vy = x_state(3);

  //pre-compute a set of terms to avoid repeated calculation
  const float c1 = (px*px) + (py*py);
  const float c2 = sqrt(c1);
  const float c3 = (c1 * c2);

  const float px_c2 = px/c2;
  const float py_c2 = py/c2; 
  const float pxvy_pyvx_c3 = (px*vy - py*vx)/c3;

  //check division by zero
  if (fabs(c1) < 0.0001) {
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    /* return a matrix initialised with 0s. */
    Hj.fill(0.0);
    return Hj;
  }

  //compute the Jacobian matrix
  Hj <<            px_c2,           py_c2,     0,     0,
                -(py/c1),         (px/c1),     0,     0,
        -py*pxvy_pyvx_c3, px*pxvy_pyvx_c3, px_c2, py_c2;

  return Hj;
}
