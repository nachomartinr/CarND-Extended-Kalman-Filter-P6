#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {

  is_initialized_ = false;

  previous_timestamp_ = 0.0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_radar_ = MatrixXd(3, 4);



  /**
  TODO:
    * Finish initializing the FusionEKF.
  */

  //Measurement model matrices. They will be updated later.
  //measurement model - laser
  H_laser_ << 1,0,0,0,
              0,1,0,0;

  //measurement model jacobian - radar
  Hj_radar_ << 0,0,0,0,
               0,0,0,0,
               0,0,0,0;


  //measurement noise covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  //measurement noise covariance matrix - radar
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  noise_ax_ = 9;
  noise_ay_ = 9;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: " << endl;

    VectorXd x = VectorXd(4);
    x << 1, 1, 1, 1;

    //initial state covariance matrix
    MatrixXd P = MatrixXd(4, 4);
    P <<  1, 0, 0,    0,
          0, 1, 0,    0,
          0, 0, 1000, 0,
          0, 0, 0,    1000;

    //state transition matrix 
    MatrixXd F = MatrixXd(4, 4);      

    //process Covariance Matrix
    MatrixXd Q = MatrixXd(4, 4);

    ekf_.Init(x, P, F, H_laser_, R_laser_, Q);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);

      float px = rho * cos(phi);
      float py = rho * sin(phi);
      // ro_dot (measurement_pack.raw_measurements_(2)) doesn't contain
      // enough information to determine vx and vy.
      float vx = 0.0;
      float vy = 0.0;

      ekf_.x_ << px, py, vx, vy;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      float px = measurement_pack.raw_measurements_(0);
      float py = measurement_pack.raw_measurements_(1);
      float vx = 0.0;
      float vy = 0.0;
      ekf_.x_ << px, py, vx, vy;
    }

    // update timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/ 1000000.0;

  /* When dt is small, skip the prediction step and effectively run multiple sequential update steps
     for the 'same' timestamp.
     This has the benefit of being a bit more efficient as it allows us to avoid the extra prediction 
     step calculations. */

  if (dt >= 0.001) {


    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    /* Update the state transition matrix */
    ekf_.F_ << 1, 0, dt, 0,
               0, 1, 0, dt,
               0, 0, 1,  0,
               0, 0, 0,  1;

    /* Update the process noise matrix */

    ekf_.Q_ <<  dt_4/4*noise_ax_, 0,               dt_3/2*noise_ax_, 0,
                0,                dt_4/4*noise_ay_, 0,                dt_3/2*noise_ay_,
                dt_3/2*noise_ax_, 0,               dt_2*noise_ax_,   0,
                0,                dt_3/2*noise_ay_, 0,                dt_2*noise_ay_;

    ekf_.Predict();
  }

  previous_timestamp_ = measurement_pack.timestamp_;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    // Update the measurement model matrix with the jacobian
    Hj_radar_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_radar_;
    /* Update the measurement noise matrix for the radar */
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
