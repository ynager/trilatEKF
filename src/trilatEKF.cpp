#include "trilatEKF.h"


TrilatEKF::TrilatEKF(VectorXd xInit, MatrixXd sensorLoc) {
    
    // initialize H
    MatrixXd H = MatrixXd(SUNFLOWER_NR, STATE_SIZE);
    
    // set state transition matrix F
    MatrixXd F = MatrixXd::Identity(STATE_SIZE, STATE_SIZE);
    
    // set measurement covariance matrix R
    MatrixXd R = MatrixXd::Identity(SUNFLOWER_NR, SUNFLOWER_NR) * VAR_Z;
    
    // set process noise covariance matrix Q
    MatrixXd Q = MatrixXd::Zero(STATE_SIZE, STATE_SIZE); //TODO remove from initialize?
    
    // set initial covariance matrix equal to noise covariance
    MatrixXd PInit = Q;
    
    // calculate initial measurement matrix H
    MatrixXd HInit = getJacobian(sensorLoc, xInit);
    
    // initialize Kalman Filter
    ekf_.initialize(xInit, PInit, F, HInit, R, Q);
    
}

TrilatEKF::~TrilatEKF() {}

void TrilatEKF::processMeasurement(const TrilatMeasurement &measurement){
    
    // ** PREDICTION **
    // get time delta dt [seconds]
    float dt = (measurement.timestamp_ - timestamp_prev_) / 1000.0;
    timestamp_prev_ = measurement.timestamp_;
    float dt2 = dt * dt;
    float dt3 = dt2 * dt;
    float dt4 = dt3 * dt;
    
    // update F matrix
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
    
    // update Q matrix
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt4/4*MEAS_NOISE_X, 0, dt3/2*MEAS_NOISE_X, 0,
               0, dt4/4*MEAS_NOISE_Y, 0, dt3/2*MEAS_NOISE_Y,
               dt3/2*MEAS_NOISE_X, 0, dt2*MEAS_NOISE_X, 0,
               0, dt3/2*MEAS_NOISE_Y, 0, dt2*MEAS_NOISE_Y;
    
    // perform prediction step
    ekf_.predict();
    
    // ** UPDATE **
    // get current jacobian
    MatrixXd H = getJacobian(measurement.sensorLoc_, ekf_.x_);
    ekf_.H_ = H;
    
    // perform update step
    ekf_.update(measurement.distance_);
    
    std::cout << "x_ = " << ekf_.x_ << std::endl;
}

MatrixXd TrilatEKF::getJacobian(const MatrixXd &sensorLoc, const VectorXd &x){
    MatrixXd H = MatrixXd::Zero(SUNFLOWER_NR, 4);
    
    for (int i = 0; i < SUNFLOWER_NR; i++) {
        double den = sqrt(pow((x(0)-sensorLoc(i,0)),2) + pow(x(1)-sensorLoc(i,1),2));
        double dh_dx = (x(0)-sensorLoc(i,1)) / den;
        double dh_dy = (x(1)-sensorLoc(i,2)) / den;
        H(i, 0) = dh_dx;
        H(i, 1) = dh_dy;
    }
    return H;
}

