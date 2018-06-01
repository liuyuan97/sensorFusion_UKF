#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.18;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.24;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  //set state dimension
  n_x_ = 5;
  //set augmented dimension
  n_aug_ = 7;
  //set spreading parameter
  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  x_ = VectorXd(n_x_);
  x_ = VectorXd::Zero(n_x_);
  //x_(2) = 5;

  P_ =  MatrixXd(n_x_, n_x_);
  P_ =  MatrixXd::Identity(n_x_, n_x_) * 0.3;

  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // initial the radar measurement variables
  Zsig_radar_ = MatrixXd(3, 2 * n_aug_ + 1);
  Z_radar_pred_ = VectorXd(3);
  S_radar_ = MatrixXd(3,3);

  H_lidar = MatrixXd(2, 5);
  R_lidar = MatrixXd(2, 2);


  H_lidar << 1.0, 0.0, 0.0, 0.0, 0.0,
             0.0, 1.0, 0.0, 0.0, 0.0;	

  R_lidar << std_laspx_ * std_laspx_, 0,
             0, std_laspy_ * std_laspy_;  

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  weights_.segment(1,2*n_aug_) =  VectorXd::Constant(2*n_aug_,1) * 0.5 /  (lambda_ + n_aug_);

  is_initialized_ = false;

  nis_file.open ("nis_data.txt");
  return;
}

UKF::~UKF() {
  nis_file.close();
  return;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if (!is_initialized_) {
   
    // first measurement
    cout << "UKF: " << endl;
  
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_(0) = meas_package.raw_measurements_(0) * cos (meas_package.raw_measurements_(1));
      x_(1) = meas_package.raw_measurements_(0) * sin (meas_package.raw_measurements_(1));

      P_(0,0) = std_radr_ * std_radr_;
      P_(1,1) = std_radr_ * std_radr_;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
      P_(0,0) = std_laspx_ * std_laspx_;
      P_(1,1) = std_laspy_ * std_laspy_;
    }

    previous_timestamp_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;

    cout << "x_ = " << x_ << endl;
    cout << "P_ = " << P_ << endl;
    return;
  }

  // return if igore radar measurements
  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && (!use_radar_)){
    return;
  }

  // return if igore lidar measurements
  if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && (!use_laser_)){
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  Prediction();
  Update(meas_package);

  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;

  return;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 */
void UKF::Prediction(void) {
  // generate the augmented sigma points
  AugmentedSigmaPoints();
  // update sigma point prediction
  SigmaPointPrediction();  
  // compute the mean and variance of the sigma points
  PredictMeanAndCovariance(); 
}


/**
 * Updates the state and the state covariance matrix using a measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::Update(MeasurementPackage meas_package) {

  double nis;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    PredictRadarMeasurement();
    nis = UpdateRadar(meas_package);
    nis_file << "R ";
  } else {
  	nis = UpdateLidar(meas_package);
  	nis_file << "L ";
  }
  nis_file << nis << endl;
  return;
}

void UKF::AugmentedSigmaPoints(void) {

  //create sigma point matrix
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd p_aug = MatrixXd(n_aug_, n_aug_);

  //create augmented mean state
  x_aug.head(n_x_)  = x_;
  x_aug.tail(n_aug_- n_x_) = VectorXd::Zero(n_aug_ - n_x_);

  //create augmented covariance matrix
  p_aug = MatrixXd::Zero(n_aug_, n_aug_);
  p_aug.block(0,0,n_x_,n_x_) = P_;
  p_aug(n_x_,n_x_) = std_a_ * std_a_;
  p_aug(n_x_+1,n_x_+1) = std_yawdd_ * std_yawdd_;

  //calculate square root of p_aug
  MatrixXd A = p_aug.llt().matrixL();
  
  float scale = sqrt(lambda_ + n_aug_);
  A = scale * A;
  
  //set sigma points as columns of matrix Xsig
  Xsig_aug_.block(0,0,n_aug_,1) = x_aug;
  Xsig_aug_.block(0,1,n_aug_,n_aug_) = A.colwise() + x_aug;
  A = -1 * A;
  Xsig_aug_.block(0,n_aug_+1,n_aug_,n_aug_) = A.colwise() + x_aug;
}

void UKF::SigmaPointPrediction(void) {

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  for (int indx = 0; indx < 2 * n_aug_ + 1; indx ++) {
      double px = Xsig_aug_(0,indx);
      double py = Xsig_aug_(1,indx);
      double v =  Xsig_aug_(2,indx);
      double psai =  Xsig_aug_(3,indx);
      double psai_d =  Xsig_aug_(4,indx);
      double a  =  Xsig_aug_(5,indx);
      double psai_dd =  Xsig_aug_(6,indx);
      if (fabs(psai_d) <= 0.001) {
          Xsig_pred_(0,indx) = Xsig_aug_(0,indx) + v * cos(psai) * delta_t + delta_t * delta_t * cos(psai) * a /2;
          Xsig_pred_(1,indx) = Xsig_aug_(1,indx) + v * sin(psai) * delta_t + delta_t * delta_t * sin(psai) * a /2;    
      } else {
          Xsig_pred_(0,indx) = Xsig_aug_(0,indx) + v /psai_d  * (sin(psai + psai_d * delta_t) - sin(psai)) + 
                              delta_t * delta_t * cos(psai) * a /2;
          Xsig_pred_(1,indx) = Xsig_aug_(1,indx) + v /psai_d  * (cos(psai) -cos(psai + psai_d * delta_t)) + 
                              delta_t * delta_t * sin(psai) * a /2;
      }
      Xsig_pred_(2,indx) = Xsig_aug_(2,indx) + a * delta_t;
      Xsig_pred_(3,indx) = Xsig_aug_(3,indx) + delta_t * psai_d + delta_t * delta_t * psai_dd / 2;
      //limit_angle(Xsig_pred_(3,indx));
      Xsig_pred_(4,indx) = Xsig_aug_(4,indx) + delta_t * psai_dd;   
  }
  return;
}


void UKF::PredictMeanAndCovariance(void) {
  //predict state mean
  x_ = VectorXd::Zero(n_x_);
  for (int indx = 0; indx < 2*n_aug_+1; indx ++) {
      x_ += weights_(indx) * Xsig_pred_.col(indx);
  }
  //predict state covariance matrix
  P_ = MatrixXd::Zero(n_x_, n_x_);
  for (int indx = 0; indx < 2*n_aug_+1; indx ++) {
      VectorXd diff = VectorXd(n_x_);
      diff = Xsig_pred_.col(indx) - x_;
      limit_angle(diff(3));
      P_ += weights_(indx)  * diff * diff.transpose();
  }
  return;
}

void UKF::PredictRadarMeasurement(void) {

  //transform sigma points into measurement space
  for (int indx = 0; indx < 2 * n_aug_ + 1; indx ++) {
      double px = Xsig_pred_(0,indx);
      double py = Xsig_pred_(1,indx);
      double v = Xsig_pred_(2,indx);
      double psai = Xsig_pred_(3,indx);
      Zsig_radar_(0,indx) = sqrt(px * px + py * py);
      Zsig_radar_(1,indx) = atan2 (py,px);
      Zsig_radar_(2,indx) = (px * cos(psai) * v + py * sin(psai) * v) / Zsig_radar_(0,indx);
  }
  //calculate mean predicted measurement
  Z_radar_pred_ = VectorXd::Zero(3);
  for (int indx = 0; indx < 2*n_aug_+1; indx ++) {
      Z_radar_pred_ += weights_(indx) * Zsig_radar_.col(indx);
  }
  
  //calculate innovation covariance matrix S
  S_radar_ = MatrixXd::Zero(3, 3);
  for (int indx = 0; indx < 2*n_aug_+1; indx ++) {
      VectorXd diff = VectorXd(3);
      diff = Zsig_radar_.col(indx) - Z_radar_pred_;
      limit_angle(diff(1));
      S_radar_ += weights_(indx)  * diff * diff.transpose();
  }
  S_radar_(0,0) += std_radr_ * std_radr_;
  S_radar_(1,1) += std_radphi_ * std_radphi_;
  S_radar_(2,2) += std_radrd_ * std_radrd_;
  return;
}

double UKF::UpdateRadar(MeasurementPackage meas_package) {
  double nis;	
  ///* the corss correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, 3);
  
  Tc = MatrixXd::Zero(n_x_, 3);
  //calculate cross correlation matrix
  for (int indx = 0; indx < 2*n_aug_+1; indx ++) {
      VectorXd x_diff = VectorXd(n_x_); 
      x_diff = Xsig_pred_.col(indx) - x_;
      limit_angle(x_diff(3));

      VectorXd z_diff = VectorXd(3); 
      z_diff = Zsig_radar_.col(indx) - Z_radar_pred_;
      limit_angle(z_diff(1));
      
      Tc += weights_(indx)  * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd S_inverse = S_radar_.inverse();
  VectorXd y = meas_package.raw_measurements_ - Z_radar_pred_;

  ///* Kalman gain
  MatrixXd K = MatrixXd(n_x_, 3);
  K = Tc * S_inverse;
  //update state mean and covariance matrix
  x_ = x_ + K * y;
  P_ = P_ - K * S_radar_ * K.transpose();

  nis = y.transpose() * S_inverse * y;
  return nis;
}

double UKF::UpdateLidar(MeasurementPackage meas_package) {
   double nis;	

   VectorXd y = VectorXd(2);

   y = meas_package.raw_measurements_ - H_lidar * x_;     

   MatrixXd Ht = H_lidar.transpose();
   MatrixXd S = H_lidar * P_ * Ht + R_lidar;
   MatrixXd S_inverse = S.inverse();
   MatrixXd K =  P_ * Ht * S_inverse;

  //new state
  x_ = x_ + (K * y);
  P_ = P_ - K * H_lidar * P_;

  nis = y.transpose() * S_inverse * y;
  return nis;
}


void UKF::limit_angle(double& angle) {
  while (angle> M_PI) angle-=2.*M_PI;
  while (angle<-M_PI) angle+=2.*M_PI;
  return;
}

