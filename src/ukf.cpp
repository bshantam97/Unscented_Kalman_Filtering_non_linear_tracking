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
  Xsig_pred_ = Eigen::MatrixXd::Zero(5, 15);

  // Initialize the weights vector
  weights_ = Eigen::VectorXd(15);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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

  // Set initialization boolean
  is_initialized_ = false;

  // Normal state dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Lambda spreading parameter
  lambda_ = 3-n_x_;

  // Measurement dimension for RADAR, measure radius, bearing and rate of change of radius
  n_z_radar_ = 3;

  n_z_lidar_ = 2;

  // Set weights
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    if (i == 0) weights_(i) = (lambda_ / (lambda_ + n_x_));
    else weights_(i) = (1 / 2*(lambda_+n_x_));
  }

  time_us_ = 0;
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
    std::cout << "Unscented Kalman Filter initialization" << std::endl;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER){

      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            0,
            0,
            0;

       // Initialize the covariance matrix. Good starting point is to initialize it with the identity
      // matrix and maintain the property of P being symmetric. Another way is to initialize it ina 
      // way that we input how much difference we expect between true state and initialzed x vector
      P_ << std_laspx_*std_laspx_,0,0,0,0,
            0,std_laspy_*std_laspy_,0,0,0,
            0,0,1,0,0,
            0,0,0,1,0,
            0,0,0,0,1;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho_ = meas_package.raw_measurements_[0];
      double phi_ = meas_package.raw_measurements_[1];
      double rho_dot_ = meas_package.raw_measurements_[2];

      // Analyze vector components
      double px_ = rho_*std::cos(phi_);
      double py_ = rho_*std::sin(phi_);

      double vx_ = rho_dot_ * std::cos(phi_);
      double vy_ = rho_dot_ * std::sin(phi_);
      double v_ = rho_dot_;   //std::sqrt(vx_*vx_ + vy_*vy_);

      x_ << px_,
            py_,
            v_,
            0,
            0;
      P_ << std_radr_*std_radr_, 0, 0, 0, 0,
            0, std_radr_*std_radr_, 0, 0, 0,
            0, 0, std_radrd_*std_radrd_, 0, 0,
            0, 0, 0, std_radphi_*std_radphi_, 0,
            0, 0, 0, 0, std_radphi_*std_radphi_;
    } else {
      std::cout << "NO MEASUREMENT RECIEVED" << std::endl;
    }

    // Input the timestamp
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  /**
   * @brief If the measurement recieved is not the first measurement
   * 
   */

  // Current - previous timestamp
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;

  // Input current time as previous time 
  time_us_ = meas_package.timestamp_;

  // Prediction step
  // Reset the values at each iteration
  if (meas_package.sensor_type_ == MeasurementPackage::LASER){

    x_ << meas_package.raw_measurements_[0],
          meas_package.raw_measurements_[1],
          0,
          0,
          0;
    P_ << std_laspx_*std_laspx_,0,0,0,0,
          0,std_laspy_*std_laspy_,0,0,0,
          0,0,1,0,0,
          0,0,0,1,0,
          0,0,0,0,1;
    Prediction(dt);
  }else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    const double rho_ = meas_package.raw_measurements_[0];
    const double phi_ = meas_package.raw_measurements_[1];
    const double rho_dot_ = meas_package.raw_measurements_[2];

    // Analyze vector components
    const double px_ = rho_*std::cos(phi_);
    const double py_ = rho_*std::sin(phi_);

    const double vx_ = rho_dot_ * std::cos(phi_);
    const double vy_ = rho_dot_ * std::sin(phi_);
    const double v_ = rho_dot_; // std::sqrt(vx_*vx_ + vy_*vy_);

    x_ << px_,
          py_,
          v_,
          0,
          0;

    P_ << std_radr_*std_radr_, 0, 0, 0, 0,
          0, std_radr_*std_radr_, 0, 0, 0,
          0, 0, std_radrd_*std_radrd_, 0, 0,
          0, 0, 0, std_radphi_*std_radphi_, 0,
          0, 0, 0, 0, std_radphi_*std_radphi_;
    Prediction(dt);
  }

  // Measurement update step
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
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
  x_aug_ << x_, 0.0, 0.0;

  // Augmented Covariance with the noise added
  Eigen::MatrixXd noise_covariance_ = Eigen::MatrixXd(2,2);
  noise_covariance_ << std_a_ * std_a_ ,0.0,
                       0.0, std_yawdd_*std_yawdd_;
  
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_.bottomRightCorner(2,2) = noise_covariance_;

  // Create square root matrix
  // This is the lower triangular matrix
  Eigen::MatrixXd A = P_aug_.llt().matrixL();

  // The first column of augmented sigma points is the mean
  Xsig_aug_.col(0) = x_aug_;

  // Generate the Sigma Points
  // Lots of indexing errors were here.Debugged
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i+1) = x_aug_ + std::sqrt(n_aug_ + lambda_)*A.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - std::sqrt(n_aug_ + lambda_)*A.col(i);
  }

  /**
   * @brief Predict the Sigma Points
   * 
   */

  for (int i = 0; i < 2*n_aug_ + 1; i++) {
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
                                                         0.0,
                                                         yaw_rate_*delta_t,
                                                         0.0;
    else process_model_ << velocity*(std::cos(yaw_)*delta_t),
                           velocity*(std::sin(yaw_)*delta_t),
                           0.0,
                           yaw_rate_*delta_t,
                           0.0;
    
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

  // State mean prediction
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    x_+= weights_(i)*Xsig_pred_.col(i);
  }

  // State Covariance prediction
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    Eigen::MatrixXd diff_ = Xsig_pred_.col(i) - x_;
    diff_(3) = std::atan2(std::sin(diff_(3)), std::cos(diff_(3)));
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
  
  // Noise matrix
  Eigen::MatrixXd additive_noise_ = Eigen::MatrixXd(n_z_lidar_,n_z_lidar_);
  additive_noise_ << std_laspx_*std_laspx_,0.0,
                     0.0, std_laspy_*std_laspy_;
  
  // Lidar measurement space. This is a matrix of predicted sigma points for the lidar
  Eigen::MatrixXd Zsig = Eigen::MatrixXd(n_z_lidar_, 2*n_aug_ + 1);

  // Measurement space conversion of Sigma points to px and py
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    Eigen::MatrixXd lidar_measurements_ = Eigen::MatrixXd(n_z_lidar_,1);
    lidar_measurements_ << Xsig_pred_.col(i)[0],
                           Xsig_pred_.col(i)[1];
    Zsig.col(i) = lidar_measurements_;
  }

  /**
   * @brief Predicted mean and covariance estimation
   * 
   */
  Eigen::VectorXd z_pred_ = Eigen::VectorXd(n_z_lidar_);

  // Measurement Covariance Matrix
  Eigen::MatrixXd S_ = Eigen::MatrixXd(n_z_lidar_,n_z_lidar_);

  // calculate mean predicted measurement
  z_pred_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
      z_pred_ += weights_(i)*Zsig.col(i);
  }
  // calculate innovation covariance matrix S
  S_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
      Eigen::MatrixXd diff_ = Zsig.col(i) - z_pred_;
      S_ += weights_(i)*(diff_)*(diff_.transpose());
  }

  // Innovation covariance for Lidar
  S_ = S_ + additive_noise_;

  /**
   * @brief Final Mean and Covariance Updates
   * 
   */
  // Create a vector of incoming LIDAR measurements
  Eigen::VectorXd z_ = Eigen::VectorXd(n_z_lidar_);
  z_ << meas_package.raw_measurements_[0],
        meas_package.raw_measurements_[1];
  
  // Cross correlation matrix calculation
  Eigen::MatrixXd Tc = Eigen::MatrixXd(n_x_, n_z_lidar_);


  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
      Eigen::MatrixXd diff_x_ = Xsig_pred_.col(i) - x_;
      Eigen::MatrixXd diff_meas_ = Zsig.col(i) - z_pred_;
      Tc += weights_(i)*diff_x_*(diff_meas_.transpose());
  }

  // Kalman Gain
  Eigen::MatrixXd K = Tc * S_.inverse();

  // Updated state
  x_ = x_ + K * (z_ - z_pred_);

  // Updated Covariance 
  P_ = P_ - K*S_*K.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  // Matrix for sigma points in measurement space
  Eigen::MatrixXd Zsig = Eigen::MatrixXd(n_z_radar_, 2*n_aug_ + 1);

  // Mean predicted measurement
  Eigen::VectorXd z_pred_ = Eigen::VectorXd(n_z_radar_);

  // Measurement covariance matrix S
  Eigen::MatrixXd S_ = Eigen::MatrixXd(n_z_radar_,n_z_radar_);

  // transform sigma points into measurement space
  Eigen::MatrixXd additive_noise_ = Eigen::MatrixXd(n_z_radar_,n_z_radar_);
  additive_noise_ << std_radr_*std_radr_ , 0.0 , 0.0,
                      0.0 , std_radphi_*std_radphi_ , 0.0,
                      0.0 , 0.0 , std_radrd_ * std_radrd_;

  // Measurent space conversion
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      Eigen::VectorXd radar_measurements_ = Eigen::VectorXd(n_z_radar_);
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
  z_pred_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
      z_pred_ += weights_(i)*Zsig.col(i);
  }
  // calculate innovation covariance matrix S
  S_.fill(0.0);
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
  Eigen::VectorXd z_ = Eigen::VectorXd(n_z_radar_);
  z_ << meas_package.raw_measurements_[0],
        meas_package.raw_measurements_[1],
        meas_package.raw_measurements_[2];

  // Calculate the cross correlation matrix
  Eigen::MatrixXd Tc = Eigen::MatrixXd(n_x_, n_z_radar_);

  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
      Eigen::MatrixXd diff_x_ = Xsig_pred_.col(i) - x_;
      // while (diff_x_(3)> M_PI) diff_x_(3)-=std::atan2(std::sin(diff_x_(3)), std::cos(diff_x_(3)));
      // while (diff_x_(3)<-M_PI) diff_x_(3)+=std::atan2(std::sin(diff_x_(3)), std::cos(diff_x_(3)));

      Eigen::MatrixXd diff_meas_ = Zsig.col(i) - z_pred_;
      // while (diff_meas_(1)> M_PI) diff_meas_(1)-=std::atan2(std::sin(diff_meas_(1)), std::cos(diff_meas_(1)));
      // while (diff_meas_(1)<-M_PI) diff_meas_(1)+=std::atan2(std::sin(diff_meas_(1)), std::cos(diff_meas_(1)));

      Tc += weights_(i)*diff_x_*(diff_meas_.transpose());
  }

  // calculate Kalman gain K;
  Eigen::MatrixXd K = Tc * S_.inverse();

  // update state mean and covariance matrix
  // z-real RADAR measurement and z_pred-predicted measurement

  x_ = x_ + K*(z_ - z_pred_);
  P_ = P_ - K*S_*K.transpose();
}