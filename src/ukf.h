#ifndef UKF_H
#define UKF_H

#include <iostream>

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  std::ofstream nis_file;	

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(void);

  /**
   * Updates the state and the state covariance matrix using a measurement
   * @param meas_package The measurement at k+1
   */
  void Update(MeasurementPackage meas_package);

private:

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* time difference between Kalman filter update
  double delta_t;

  ///* previous timestamp
  long previous_timestamp_;

  ///* aumentation sigma points
  MatrixXd Xsig_aug_;

  ///* radar sigma points in measurement space
  MatrixXd Zsig_radar_;

  ///* radar mean predicted measurement
  VectorXd Z_radar_pred_;
  
  ///* radar measurement covariance matrix S
  MatrixXd S_radar_;

  MatrixXd H_lidar;
  MatrixXd R_lidar;

  void AugmentedSigmaPoints(void);
  void SigmaPointPrediction(void);
  void PredictMeanAndCovariance(void);
  void PredictRadarMeasurement(void);
  double UpdateRadar (MeasurementPackage meas_package);
  double UpdateLidar (MeasurementPackage meas_package);
  void limit_angle(double& angle);
};

#endif /* UKF_H */
