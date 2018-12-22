#include "trilatEKF.h"


TrilatEKF::TrilatEKF(VectorXd xInit, MatrixXd sensorLoc) {
    
    MatrixXd H = MatrixXd(SUNFLOWER_NR, STATE_SIZE);                            // initialize H
    MatrixXd F = MatrixXd::Identity(STATE_SIZE, STATE_SIZE);                    // set state transition matrix F
    MatrixXd R = MatrixXd::Identity(SUNFLOWER_NR, SUNFLOWER_NR) * VAR_Z;        // set measurement covariance matrix R
    MatrixXd Q = MatrixXd::Zero(STATE_SIZE, STATE_SIZE);                        // TODO remove from initialize?
    
    MatrixXd PInit = Q;                                 // set initial covariance matrix equal to noise covariance
    std::cout << "xInit: " << xInit << std::endl;
    ekf_.x_ = xInit; // TODO HACK
    MatrixXd HInit = getJacobian(sensorLoc, xInit);     // calculate initial measurement matrix H
    
    // initialize Kalman Filter
    ekf_.initialize(xInit, PInit, F, HInit, R, Q);
    
}

TrilatEKF::~TrilatEKF() {}

void TrilatEKF::processMeasurement(const TrilatMeasurement &measurement){
    
    // ** PREDICTION **
    float dt = (measurement.timestamp_ - timestamp_prev_) / 1000.0; // get time delta dt [seconds]
    float dt2 = dt * dt;
    float dt3 = dt2 * dt;
    float dt4 = dt3 * dt;
    timestamp_prev_ = measurement.timestamp_;
    
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
    MatrixXd H = getJacobian(measurement.sensorLocs_, ekf_.x_);
    ekf_.H_ = H;
    
    // perform update step
    ekf_.update(measurement.distances_);
}

MatrixXd TrilatEKF::getJacobian(const MatrixXd &sensorLoc, const VectorXd &x){
    MatrixXd sensorLoc2 = MatrixXd(3,2); // initial sensor location
    sensorLoc2 << 0.0, 0.0,
    5.0, 0.0,
    0.0, 5.0;
    
    MatrixXd H = MatrixXd::Zero(SUNFLOWER_NR, 4);
    
    for (int i = 0; i < SUNFLOWER_NR; i++) {
        // std::cout << "x: " << x(0) << ", " << x(1) << std::endl;
        // std::cout << "sensor loc: " << sensorLoc2 << std::endl;
        //std::cout << "sensor loc: " << sensorLoc(i,0) << ", " << sensorLoc(i,1) << std::endl;
        double den = sqrt(pow((x(0)-sensorLoc2(i,0)),2) + pow(x(1)-sensorLoc2(i,1),2));
        if (den < 0.00001) {
            std::cout << "Error: Division by Zero in getJacobian()" << std::endl;
            return H;
        }
        double dh_dx = (x(0)-sensorLoc2(i,0)) / den;
        double dh_dy = (x(1)-sensorLoc2(i,1)) / den;
        H(i, 0) = dh_dx;
        H(i, 1) = dh_dy;
    }
    return H;
}

