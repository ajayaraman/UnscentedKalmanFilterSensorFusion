#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  n_x_ = 5;
  n_aug_ = 7;


  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = pow(M_PI / 8.0 , 2);

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  H_laser_ = MatrixXd(2, n_x_);
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;  

  R_laser_ = MatrixXd(2,2);
  //measurement covariance matrix - laser
  R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.);

  weights_ = VectorXd(2*n_aug_+1);

  //set weights
  weights_[0] = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2*n_aug_+1; i++)
  {
    weights_[i] = 1 / (2* (lambda_ + n_aug_)); 
  }                 
           
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

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

    VectorXd meas = meas_package.raw_measurements_;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //rho * cos(phi), rho * sin(phi)
      x_(0) = meas(0) * cos(meas(1));
      x_(1) = meas(0) * sin(meas(1));
      x_(3) = meas(2); //We also measure phi which is in the state vector
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      // x, y 
      x_(0) = meas(0);
      x_(0) = meas(1);
    }

    //setup previous timestamp to be the timestamp of the first packet so that dt is initialized correctly
    time_us_ = meas_package.timestamp_;
   
    // done initializing, no need to predict or update
    is_initialized_ = true;

    //std::cout << "Is initialized" << std::endl;
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
    //Debug
    //std::cout << "Update after Laser x = " << x_ << std::endl;
    //std::cout << "Update after Laser P = " << P_ << std::endl;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
    //Debug
    //std::cout << "Update after Radar x = " << x_ << std::endl;
    //std::cout << "Update after Radar P = " << P_ << std::endl;    
  }

  time_us_ = meas_package.timestamp_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  
  //Compute Augmented state
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug[5] = 0;
  x_aug[6] = 0 ;
  
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner( n_x_ , n_x_ ) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;
  
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  
  //create augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //set first column of sigma point matrix
  Xsig_aug.col(0)  = x_aug;

  //set remaining sigma points
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)     = x_aug + sqrt(lambda_+n_aug_) * A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * A.col(i);
  }

 //predict sigma points according to the motion model CTRV
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //Predicted Mean and Covariance

  //predict state mean
  x_.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
      x_ += weights_[i] * Xsig_pred_.col(i);
  }
  
  //predict state covariance matrix
  P_.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
      VectorXd xdiff = Xsig_pred_.col(i) - x_;
      while(xdiff(3) > M_PI ) xdiff(3) -= 2 * M_PI;
      while(xdiff(3) < -M_PI) xdiff(3) += 2 * M_PI;
      P_ += weights_[i] * xdiff * xdiff.transpose();
  }

  //Debug
  //std::cout << "x = " << x_ << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  MatrixXd Ht = H_laser_.transpose();
  VectorXd z_diff = meas_package.raw_measurements_ - H_laser_ * x_;
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd innovation = S.inverse();
  MatrixXd K = P_ * Ht * innovation;
		
  x_ = x_ + K * z_diff;

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;

  //Calculate NIS Score
  double NIS = z_diff.transpose() * innovation * z_diff;
  std::cout << "NIS Radar = " << NIS << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //create matrix for sigma points in measurement space
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  //transform sigma points into measurement space
  for(int i = 0 ; i < 2 * n_aug_ + 1; i ++)
  {
      double px = Xsig_pred_(0, i);
      double py = Xsig_pred_(1, i);
      double v =  Xsig_pred_(2, i);
      double psi = Xsig_pred_(3, i);
      double psidot = Xsig_pred_(4, i);
      
      Zsig(0, i) = sqrt(px * px + py * py);
      Zsig(1, i) = atan2(py, px); 
      Zsig(2, i) = (px * cos(psi) + py * sin(psi)) * v / (Zsig(0,i) + 0.00001);
  }
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1 ; i++)
  {
      z_pred += weights_[i] * Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  
  S.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++)
  {
      VectorXd zdiff = Zsig.col(i) - z_pred;
      S += weights_[i] * (zdiff * zdiff.transpose());
  }
  
  MatrixXd R(n_z, n_z);
  
  R.fill(0);
  //Diagonal matrix with measurement noise
  R(0,0) = std_radr_ * std_radr_;
  R(1,1) = std_radphi_ * std_radphi_;
  R(2,2) = std_radrd_ * std_radrd_;
  S += R; 

  //Compute Kalman gain and update states and covariance matrix

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0);
  //calculate cross correlation matrix
  for(int i = 0; i < 2 * n_aug_ + 1; i++)
  {
      VectorXd xdiff = Xsig_pred_.col(i) - x_;
      while(xdiff(3) > M_PI) xdiff(3) -= 2 * M_PI;
      while(xdiff(3) < -M_PI) xdiff(3) += 2 * M_PI;
      VectorXd zdiff = Zsig.col(i) - z_pred;
      while(zdiff(1) > M_PI) zdiff(1) -= 2 * M_PI;
      while(zdiff(1) < -M_PI) zdiff(1) += 2 * M_PI;
      Tc += weights_[i] * xdiff * zdiff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd innovation = S.inverse();

  MatrixXd K = Tc * innovation;
  //update state mean and covariance matrix
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  //Calculate NIS Score
  double NIS = z_diff.transpose() * innovation * z_diff;
  std::cout << "NIS Radar = " << NIS << std::endl;
}
