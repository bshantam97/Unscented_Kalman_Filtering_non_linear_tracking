#include "ukf.h"
#include "Eigen/Dense"

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = Eigen::VectorXd(5);

  // initial covariance matrix
  P_ = Eigen::MatrixXd(5, 5);

  // Predicted Sigma points matrix
  Eigen::MatrixXd Xsig_pred_ = Eigen::MatrixXd(n_x_, 2*n_aug_ + 1);

  // Initialize the weights vector
  Eigen::VectorXd weights_ = Eigen::VectorXd(2*n_aug_ + 1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
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

  // Normal state dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Lambda spreading parameter
  lambda_ = 3 - n_aug_;

  // Measurement dimension for RADAR, measure radius, bearing and rate of change of radius
  n_z_ = 3;

  // Initialize the covariance matrix. Good starting point is to initialize it with the identity
  // matrix and maintain the property of P being symmetric. Another way is to initialize it ina 
  // way that we input how much difference we expect between true state and initialzed x vector

  P_ << 1,0,0,0,0,
        0,1,0,0,0,
        0,0,1,0,0,
        0,0,0,1,0,
        0,0,0,0,1;
  
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  // Initialize the px, py and can play around with the other hyperparameters to see what works best
  // RADAR velocity is not the same as CTRV velocity as RADAR velocity is measured from 
  // the vehicle's perspective and the CTRV velocity is from the objects perspective
  // Radar velocity is measured from the autonomous vehicle's perspective,In the CTRV model, the velocity is from the object's perspective
  // If you drew a straight line from the vehicle to the bicycle, radar measures the velocity along that line
  // which in this case is the bicycle; the CTRV velocity is tangential to the circle 
  // along which the bicycle travels. Therefore, you cannot directly use the 
  // radar velocity measurement to initialize the state vector.
  if (!is_initialized_) {
    x_ << meas_package.raw_measurements_[0],
          meas_package.raw_measurements_[1],
          0, // Velocity magnitude
          0, // yaw
          0; // yaw rate
    is_initialized_ = true;
  } 
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  /**
   * @brief Sigma Point Generation
   * 
   */

  // Augmented mean vector
  Eigen::VectorXd x_aug_ = Eigen::VectorXd(7);

  // Augmented State covariance
  Eigen::MatrixXd P_aug_ = Eigen::MatrixXd(7,7);

  // Generated Sigma Point Matrix
  Eigen::MatrixXd Xsig_aug_ = Eigen::MatrixXd(n_aug_, 2*n_aug_ + 1);

  // Create the augmented mean state
  x_aug_ << x_, 0, 0;

  // Augmented Covariance with the noise added
  Eigen::MatrixXd noise_covariance_ = Eigen::MatrixXd(2,2);
  noise_covariance_ << std_a_ * std_a_ ,0,
                       0, std_yawdd_*std_yawdd_;
  
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_.bottomRightCorner(2,2) = noise_covariance_;

  // Create square root matrix
  // This is the lower triangular matrix
  Eigen::MatrixXd A = P_aug_.llt().matrixL();

  // The first column of augmented sigma points is the mean
  Xsig_aug_.col(0) = x_aug_;

  // Generate the Sigma Points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    Xsig_aug_.col(i) = x_aug_ + std::sqrt(n_aug_ + lambda_)*A.col(i-1);
    Xsig_aug_.col(i+n_aug_) = x_aug_ - std::sqrt(n_aug_ + lambda_)*A.col(i-1);
  }

  /**
   * @brief Predict the Sigma Points
   * 
   */

  for (int i = 0; i < Xsig_aug_.cols(); i++) {
    // Coefficients
    double velocity = Xsig_aug_.col(i)[2];
    double yaw_ = Xsig_aug_.col(i)[3];
    double yaw_rate_ = Xsig_aug_.col(i)[4];
    double acc_noise_ = Xsig_aug_.col(i)[5];
    double yaw_acc_noise_ = Xsig_aug_.col(i)[6];

    Eigen::MatrixXd process_model_ = Eigen::MatrixXd(5,1);
    // Avoid division by zero
    if (std::fabs(yaw_rate_) > 0.001) process_model_ << (velocity / yaw_rate_)*(std::sin(yaw_ + yaw_rate_*delta_t) - std::sin(yaw_)),
                                                        (velocity / yaw_rate_)*(-std::cos(yaw_ + yaw_rate_*delta_t) + std::cos(yaw_)),
                                                        0,
                                                        yaw_rate_*delta_t,
                                                        0;
    else process_model_ << velocity*(std::cos(yaw_)*delta_t),
                           velocity*(std::sin(yaw_)*delta_t),
                           0,
                           0,
                           0;
    
    Eigen::MatrixXd process_noise_ = Eigen::MatrixXd(5,1);
    process_noise_ << 0.5*(delta_t*delta_t)*std::cos(yaw_)*acc_noise_,
                      0.5*(delta_t*delta_t)*std::sin(yaw_)*acc_noise_,
                      delta_t*acc_noise_,
                      0.5*(delta_t*delta_t)*yaw_acc_noise_,
                      delta_t*yaw_acc_noise_;
    // Write predicted sigma points into the right column
    Xsig_pred_.col(i) = Xsig_aug_.col(i).head(n_x_) + process_model_ + process_noise_;
  }

  /**
   * @brief Predict the mean and covariance
   * 
   */

  // Set weights
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    if (i == 0) weights_(i) = (lambda_ / (lambda_ + n_aug_));
    else weights_(i) = (0.5 / (lambda_+n_aug_));
  }

  // State mean prediction
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    x_+= weights_(i)*Xsig_pred_.col(i);
  }

  // State Covariance prediction
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    Eigen::MatrixXd diff_ = Xsig_pred_.col(i) - x_;
    while (diff_(3) > M_PI) diff_(3) -= 2.*M_PI;
    while (diff_(3) < -M_PI) diff_(3) += 2.*M_PI;
    P_ += weights_(i)*(diff_)*(diff_.transpose());
  } 
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  // Matrix for sigma points in measurement space
  Eigen::MatrixXd Zsig = Eigen::MatrixXd(n_z_, 2*n_aug_ + 1);

  // Mean predicted measurement
  Eigen::VectorXd z_pred_ = Eigen::VectorXd(n_z_);

  // Measurement covariance matrix S
  Eigen::MatrixXd S_ = Eigen::MatrixXd(n_z_,n_z_);

  // transform sigma points into measurement space
  Eigen::MatrixXd additive_noise_ = Eigen::MatrixXd(3,3);
  additive_noise_ << std_radr_*std_radr_ , 0 , 0,
                      0 , std_radphi_*std_radphi_ , 0,
                      0 , 0 , std_radrd_ * std_radrd_;

  // Measurent space conversion
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      Eigen::MatrixXd radar_measurements_ = Eigen::MatrixXd(3,1);
      double px_ = Xsig_pred_.col(i)[0];
      double py_ = Xsig_pred_.col(i)[1];
      double yaw_ = Xsig_pred_.col(i)[3];
      double velocity_ = Xsig_pred_.col(i)[2];
      
      double rho_ = std::sqrt(std::pow(px_, 2) + std::pow(py_, 2));
      double phi_ = std::atan2(py_, px_);
      double rho_dot_ = (px_*std::cos(yaw_)*velocity_ + py_*std::sin(yaw_)*velocity_) / std::sqrt(std::pow(px_,2) + std::pow(py_, 2));

      radar_measurements_ <<  rho_,
                              phi_,
                              rho_dot_;
      Zsig.col(i) = radar_measurements_;
    }
  // calculate mean predicted measurement
  for (int i = 0; i < 2*n_aug_+1; i++) {
      z_pred_ += weights_(i)*Zsig.col(i);
  }
  // calculate innovation covariance matrix S
  for (int i = 0; i < 2*n_aug_+1; i++) {
      Eigen::MatrixXd diff_ = Zsig.col(i) - z_pred_;
      S_ += weights_(i)*(diff_)*(diff_.transpose());
  }

  // Innovation Covariance for Radar
  S_ = S_ + additive_noise_;  

  /**
   * @brief State Matrix and Covariance Matrix Update
   * 
   */
  // Now use the incoming RADAR measurements to update the state 
  Eigen::VectorXd z_ = Eigen::VectorXd(n_z_);
  z_ << meas_package.raw_measurements_[0],
        meas_package.raw_measurements_[1],
        meas_package.raw_measurements_[2];

  // Calculate the cross correlation matrix
  Eigen::MatrixXd Tc = Eigen::MatrixXd(n_x_, n_z_);

  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
      Eigen::MatrixXd diff_x_ = Xsig_pred_.col(i) - x_;
      Eigen::MatrixXd diff_meas_ = Zsig.col(i) - z_pred_;
      Tc += weights_(i)*diff_x_*(diff_meas_.transpose());
  }

  // calculate Kalman gain K;
  Eigen::MatrixXd K = Tc * S_.inverse();

  // update state mean and covariance matrix
  // z-real RADAR measurement and z_pred-predicted measurement

  x_ = x_ + K*(z_ - z_pred_);
  P_ = P_ - K*S_*K.transpose();
}