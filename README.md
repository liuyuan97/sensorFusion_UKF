# Unscented Kalman Filter Project
Self-Driving Car Engineer Nanodegree Program

In this project, an Unscented Kalman Filter is utilized to estimate the state of a moving object of interest with noisy lidar and radar measurements. Passing the project requires obtaining RMSE values that are lower that the tolerance outlined in the project rubric. 

[//]: # (Image References)

[image1]: ./pics/lidar_nis.png "lidar Image"

[image2]: ./pics/radar_nis.png "radar Image"


## Requirements

* The Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases).  
* Set up and install [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for either Linux or Mac systems

## Build & Run

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./UnscentedKF

## Components

* main.cpp -  reads in data, calls a function to run the Unscented Kalman filter, calls a function to calculate RMSE
* ukf.cpp - initializes the Unscented Kalman filter, calls the predict and update function, defines the predict and update functions
* tools.cpp- function to calculate RMSE

## Process Noise Parameter Tuning

For the CTRV model, two parameters define the process noise:  
* longitudinal acceleration noise  
* yaw acceleration noise  

These two parameters have to choose properly to achieve the optimum tracking of the UKF.  A good to measure the noise level is to compare Normalized Innovation Square (NIS) from the sensor measurement against a threshold which is determined by the measurement dimension.  Most of NIS value should be under that threshold.  The following figure gives the NIS of the Lidar and Radar sensors.  It shows that the process noise level are close the actual values.

![alt text][image1]
![alt text][image2]