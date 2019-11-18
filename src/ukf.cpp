#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  is_initialized_ = false;

  time_us_ = 0;
  
  n_x_ = 5;

  n_aug_ = 7;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  /**
   * Initialization
   */
  
  if(!is_initialized_) {

    // first measurement
    x_ << 1, 1, 1, 1, 0.1;

    // init covariance matrix
    P_ << 0.15,    0, 0, 0, 0,
          0, 0.15, 0, 0, 0,
          0,    0, 1, 0, 0,
          0,    0, 0, 1, 0,
          0,    0, 0, 0, 1;

    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double ro = meas_package.raw_measurements_(0);
      double theta = meas_package.raw_measurements_(1);
      double ro_dot = meas_package.raw_measurements_(2);
      x_(0) = ro * cos(theta);
      x_(1) = ro * sin(theta);
    }

    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double px = meas_package.raw_measurements_(0);
      double py = meas_package.raw_measurements_(1);
      x_(0) = px;
      x_(1) = py;
    }

  time_us_ = meas_package.timestamp_;
  is_initialized_ = true;
  return;
  }
  
  /**
   * Prediction
   */
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);
  
  if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
  
    UpdateLidar(meas_package);

  }

  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
  
    UpdateRadar(meas_package);
  
  }
  
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  
  // 1. Generating Augmented Sigma Points
  lambda_ = 3 - n_aug_;

  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(7, 7);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  MatrixXd Q = MatrixXd(2,2);
  Q << std_a_ * std_a_, 0.0, 
       0.0, std_yawdd_ * std_yawdd_;

  x_aug.head(5) = x_;
  // louis modify start
  x_aug(5) = 0;
  x_aug(6) = 0;
  P_aug.fill(0.0);
  // louis modify end
  P_aug.topLeftCorner(5,5) = P_;
  P_aug.bottomRightCorner(2,2) = Q;

  MatrixXd B = P_aug.llt().matrixL();
  MatrixXd Sig3 = MatrixXd::Zero(n_aug_,n_aug_);
  MatrixXd Sig4 = MatrixXd::Zero(n_aug_,n_aug_);  

  for (int i = 0; i < B.cols(); i++) {
    Sig3.col(i) = sqrt(n_aug_ + lambda_) * B.col(i) + x_aug;
    Sig4.col(i) = x_aug - sqrt(n_aug_ + lambda_) * B.col(i);
  }

  Xsig_aug.col(0) = x_aug;

  for(int i = 1; i < n_aug_ + 1; i++) {
    Xsig_aug.col(i) = Sig3.col(i-1);
    Xsig_aug.col(i+7) = Sig4.col(i-1);
  }

  // std::cout << Xsig_aug << std::endl;
  // std::cout << std::endl;

  // 2. Predict Sigma Points  
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
        py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
        px_p = p_x + v * delta_t * cos(yaw);
        py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  // 3. Predict Mean And Covariance
  // create vector for weights
  VectorXd weights = VectorXd(2 * n_aug_ + 1);

  // create vector for predicted state
  VectorXd x_p = VectorXd(n_x_);

  // create covariance matrix for prediction
  MatrixXd P_p = MatrixXd(n_x_, n_x_);  

  // set weights
  weights(0) = lambda_ / (lambda_ + n_aug_);  

  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights(i) = weight;
  }  
  
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     x_ = x_ + weights(i) * Xsig_pred_.col(i);
  }

  P_.fill(0.0);
  // predict state covariance matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights(i) * x_diff * x_diff.transpose() ;    
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  int n_z = 2;

  // set vector for weights
  VectorXd weights = VectorXd(2 * n_aug_ + 1);  

  // set weights
  weights(0) = lambda_ / (lambda_ + n_aug_);  

  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights(i) = weight;
  }  

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    // measurement model
    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;
  }
  
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights(i) * Zsig.col(i);
  }
  
  S.fill(0.0);
  // calculate innovation covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + z_diff * z_diff.transpose() * weights(i);
  }

  // add measurement noise
  S = S + R;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_,n_z);
  MatrixXd K = MatrixXd(n_x_,n_z);

  Tc.fill(0.0);
  // calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;
    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  K = Tc * S.inverse();
  
  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;

  P_ = P_ - K * S * K.transpose();

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  
  // 1. Predict Radar Measurement
  
  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3; 

  // set vector for weights
  VectorXd weights = VectorXd(2 * n_aug_ + 1);  
  // set weights
  weights(0) = lambda_ / (lambda_ + n_aug_);  

  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights(i) = weight;
  }  

  
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_radr_ * std_radr_, 0, 0,
       0, std_radphi_ * std_radphi_, 0,
       0, 0, std_radrd_ * std_radrd_;

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double yaw_dot = Xsig_pred_(4,i);

    Zsig(0,i) = sqrt(px * px + py * py);
    Zsig(1,i) = atan2(py , px);
    Zsig(2,i) = (px * cos(yaw) * v + py * sin(yaw) * v) / sqrt(px * px + py* py);
  }
  
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  // calculate innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // anlge normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
    S = S + z_diff * z_diff.transpose() * weights(i);
  }

  // add measurement noise
  S = S + R;

  // 2. Update State

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_,n_z);
  MatrixXd K = MatrixXd(n_x_,n_z);
  
  Tc.fill(0.0);
  // calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(1) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(1) += 2.*M_PI;
    
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
    
    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  K = Tc * S.inverse();

  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;


  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;

  P_ = P_ - K * S * K.transpose();
}
