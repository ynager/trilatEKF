#include "trilatEKF.h"


TrilatEKF::TrilatEKF(VectorXd xInit, MatrixXd sensorLoc) {
    
    timestamp_prev_ = 0;
    MatrixXd F = MatrixXd::Identity(STATE_SIZE, STATE_SIZE);                    // set state transition matrix F
    MatrixXd R = MatrixXd::Identity(SUNFLOWER_NR, SUNFLOWER_NR) * VAR_Z;        // set measurement covariance matrix R
    MatrixXd Q = MatrixXd::Zero(STATE_SIZE, STATE_SIZE);                        // TODO remove from initialize?
    
    MatrixXd PInit = MatrixXd::Identity(STATE_SIZE, STATE_SIZE)*0.1;  // initial state is known
    std::cout << "xInit: " << xInit << std::endl;
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
    float dt5 = dt4 * dt;
    timestamp_prev_ = measurement.timestamp_;
    
    // update F matrix
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
    
    // update Q matrix
    ekf_.Q_ = MatrixXd(4, 4);
    //        ekf_.Q_ << dt4/4*MEAS_NOISE_X, 0, dt3/2*MEAS_NOISE_X, 0,
    //                   0, dt4/4*MEAS_NOISE_Y, 0, dt3/2*MEAS_NOISE_Y,
    //                   dt3/2*MEAS_NOISE_X, 0, dt2*MEAS_NOISE_X, 0,
    //                   0, dt3/2*MEAS_NOISE_Y, 0, dt2*MEAS_NOISE_Y;
    
    ekf_.Q_ <<  dt5/20*MEAS_NOISE_X, 0, dt4/8*MEAS_NOISE_X, 0,
    0, dt5/20*MEAS_NOISE_Y, 0, dt4/8*MEAS_NOISE_Y,
    dt4/8*MEAS_NOISE_X, 0, dt3/3*MEAS_NOISE_X, 0,
    0, dt4/8*MEAS_NOISE_Y, 0, dt3/3*MEAS_NOISE_Y;
    //    ekf_.Q_ = ekf_.Q_ * 0.00006;
    
    
    // perform prediction step
    ekf_.predict();
    
    // ** UPDATE **
    // get current jacobian
    MatrixXd H = getJacobian(measurement.sensorLocs_, ekf_.x_);
    // TODO set H row to zero if sensor out of range
    ekf_.H_ = H;
    
    // perform update step
    ekf_.update(measurement.distances_);
}

MatrixXd TrilatEKF::getJacobian(const MatrixXd &sensorLoc, const VectorXd &x){
    
    MatrixXd H = MatrixXd::Zero(SUNFLOWER_NR, 4);
    for (int i = 0; i < SUNFLOWER_NR; i++) {
        //std::cout << "x: " << x(0) << ", " << x(1) << std::endl;
        //std::cout << "sensor loc: " << sensorLoc << std::endl;
        //std::cout << "sensor loc: " << sensorLoc(i,0) << ", " << sensorLoc(i,1) << std::endl;
        double den = sqrt(pow((x(0)-sensorLoc(i,0)),2) + pow((x(1)-sensorLoc(i,1)),2));
        if (den < 0.001) {
            std::cout << "Error: Division by zero in getJacobian()\n";
            return H;
        }
        double dh_dx = (x(0)-sensorLoc(i,0)) / den;
        double dh_dy = (x(1)-sensorLoc(i,1)) / den;
        H(i, 0) = dh_dx;
        H(i, 1) = dh_dy;
        
    }
    return H;
}

TrilatMeasurement TrilatEKF::matchMeasurements(std::vector<Measurement> mvec){
    TrilatMeasurement tm;
    tm.timestamp_ = mvec[0].timestamp_;
    tm.sensorLocs_ = MatrixXd(3,2); // TODO magic var
    tm.sensorLocs_ << mvec[0].sensorLoc_.transpose(), mvec[1].sensorLoc_.transpose(), mvec[2].sensorLoc_.transpose();
    tm.distances_ = Vector3d();
    tm.distances_ << mvec[0].distance_, mvec[1].distance_, mvec[2].distance_;
    //std::cout << "locs: " << tm.sensorLocs_ << std::endl;
    //std::cout << "dist: " << tm.distances_ << std::endl;
    return tm;
}

