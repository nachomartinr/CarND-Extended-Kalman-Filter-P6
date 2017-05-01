#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#define CHOLESKY_DECOMP_INVERSE

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
  long x_size = x_in.size();
  I_ = MatrixXd::Identity(x_size, x_size);
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

/* Get the inverse of Matrix S */
static inline MatrixXd Inverse(const MatrixXd & S) {

#if defined(CHOLESKY_DECOMP_INVERSE)
  /* Assuming that S is symetric positive definite, we can apply the 
  Cholesky decomposition as a more efficient way of calculating the inverse.
  In Eigen, this is done with Si = S.llt().solve(I). 
  It might fail if numerical errors cause the S matrix to
  lose its positive definite properties. */
  return S.llt().solve(MatrixXd::Identity(S.rows(), S.cols()));

#else
  return S.inverse();
#endif

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd PHt = P_ * H_.transpose(); /* reuse P * H' */
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Si = Inverse(S);
  MatrixXd K = PHt * Si;
  
  x_ = x_ + K * y;  
  P_ = (I_ - K * H_) * P_;
}

/* Constraint angle to [-PI, PI] */
static inline float WrapAngle(const float & a) {

  float wrapped;
  if (( a >= -M_PI ) && (a <= M_PI )) {
    wrapped = a;
  } 
  else if ( a < 0) {
    wrapped = fmod(a - M_PI, 2*M_PI) + M_PI;
  }
  else {
    wrapped = fmod(a + M_PI, 2*M_PI) - M_PI;
  }

  return wrapped;
}


void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  // Apply h(x) to calculate predicted measurement z_pred = (rho, phi, rho_dot) from the predicted state
  // x_ = (px, py, vx, vy)
  float rho = sqrt(x_[0]*x_[0] + x_[1]*x_[1]);
  float phi, rho_dot;

  if (fabs(rho > 0.0001)) {
    phi = atan2(x_[1], x_[0]);
    rho_dot = ((x_[0]*x_[2] + x_[1]*x_[3]) / rho);
  } 
  else {
    // px and py are very small. Make phi and rho_dot 0.
    phi = 0.0;
    rho_dot = 0.0;
  }
  
  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;

  VectorXd y = z - z_pred;
  /* Constraint angle to [-PI, PI] */
  y(1) = WrapAngle(y(1));

  MatrixXd PHt = P_ * H_.transpose(); /* reuse P * H' */
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Si = Inverse(S);  
  MatrixXd K = PHt * Si;
  
  x_ = x_ + K * y;  
  P_ = (I_ - K * H_) * P_;

}
