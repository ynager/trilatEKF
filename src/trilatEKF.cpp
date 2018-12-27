#include "trilatEKF.h"


TrilatEKF::TrilatEKF(VectorXd xInit, MatrixXd sensorLoc) {
    
    // TODO improve initialization
    timestamp_prev_ = 0;
    sensorLoc_ = sensorLoc;
    MatrixXd F = MatrixXd::Identity(STATE_SIZE, STATE_SIZE);                    // set state transition matrix F
    MatrixXd R = MatrixXd::Identity(SUNFLOWER_NR, SUNFLOWER_NR) * VAR_Z;        // set measurement covariance matrix R
    MatrixXd Q = MatrixXd::Identity(STATE_SIZE, STATE_SIZE);                    // set process noise covariance
    
    MatrixXd PInit = MatrixXd::Zero(STATE_SIZE, STATE_SIZE);                    // only position of initial state is known
    PInit(2,2) = 10;
    PInit(3,3) = 10;
    std::cout << "xInit: " << xInit.transpose() << std::endl;
    MatrixXd HInit = getJacobian(sensorLoc, xInit);                             // calculate initial measurement matrix H
    
    ekf_.initialize(xInit, PInit, F, HInit, R, Q);                              // initialize Kalman Filter with xInit
    
}

TrilatEKF::~TrilatEKF() {}

void TrilatEKF::processMeasurements(std::vector<TrilatMeasurement>* measurements){
    
    // ** PREDICTION **
    double dt = ((*measurements)[0].timestamp_ - timestamp_prev_) / 1000.0; // get time delta dt [seconds]
    double dt2 = dt * dt;
    double dt3 = dt2 * dt;
    double dt4 = dt3 * dt;
    double dt5 = dt4 * dt;
    timestamp_prev_ = (*measurements)[0].timestamp_;
    
    // update F matrix
    switch(STATE_SIZE) {
            
        case 6:
            ekf_.F_ <<  1, 0, dt, 0, dt2/2, 0,
            0, 1, 0, dt, 0, dt2/2,
            0, 0, 1, 0, dt, 0,
            0, 0, 0, 1, 0, dt,
            0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 1;
        case 4:
            ekf_.F_(0, 2) = dt;
            ekf_.F_(1, 3) = dt;
            break;
            
    }
    
    // update Q matrix
    ekf_.Q_ = MatrixXd(STATE_SIZE, STATE_SIZE);
    switch(STATE_SIZE) {
            
        case 6: // PVA
            ekf_.Q_ <<  dt5/20*MEAS_NOISE_X, 0, dt4/8*MEAS_NOISE_X, 0, dt3/6*MEAS_NOISE_X, 0,
            0, dt5/20*MEAS_NOISE_Y, 0, dt4/8*MEAS_NOISE_Y, 0, dt3/6*MEAS_NOISE_Y,
            dt4/8*MEAS_NOISE_X, 0, dt3/3*MEAS_NOISE_X, 0, dt2/2*MEAS_NOISE_X, 0,
            0, dt4/8*MEAS_NOISE_Y, 0, dt3/3*MEAS_NOISE_Y, 0, dt2/2*MEAS_NOISE_X,
            dt3/3*MEAS_NOISE_X, 0, dt2/2*MEAS_NOISE_X, 0, dt*MEAS_NOISE_X, 0 ,
            0, dt3/3*MEAS_NOISE_Y, 0, dt2/2*MEAS_NOISE_Y, 0, dt*MEAS_NOISE_Y;
            break;
            
        case 4: // PV
            ekf_.Q_ << dt4/4*MEAS_NOISE_X, 0, dt3/2*MEAS_NOISE_X, 0,
            0, dt4/4*MEAS_NOISE_Y, 0, dt3/2*MEAS_NOISE_Y,
            dt3/2*MEAS_NOISE_X, 0, dt2*MEAS_NOISE_X, 0,
            0, dt3/2*MEAS_NOISE_Y, 0, dt2*MEAS_NOISE_Y;
            break;
    }
    
    
    // perform prediction step
    // std::cout << "x: " << ekf_.x_(0) << ", " << ekf_.x_(1) << std::endl;
    ekf_.predict();
    // std::cout << "x predicted: " << ekf_.x_(0) << ", " << ekf_.x_(1) << std::endl;
    
    // ** UPDATE **
    MatrixXd H = getJacobian((*measurements)[0].sensorLocs_, ekf_.x_);     // get current jacobian
    ekf_.H_ = H;
    // std::cout << "H: " << H << "\n\n";
    // TODO set H row to zero if sensor out of range
    
    // put distances into Matrix
    MatrixXd distanceMatrix(SUNFLOWER_NR, (*measurements).size());
    for (int i = 0; i < (*measurements).size(); i++)
        distanceMatrix.col(i) = (*measurements)[i].distances_;
    
    // perform update step
    int idx = ekf_.updateMahalanobis(distanceMatrix);
    
    // remove measurements from vector containing distance readings that were selected for this object
    TrilatMeasurement a = (*measurements)[(idx + 4) % 8]; // measurements for other object known since only two readings per sensor
    (*measurements).clear();
    (*measurements).push_back(a);
}

MatrixXd TrilatEKF::getJacobian(const MatrixXd &sensorLoc, const VectorXd &x){
    
    MatrixXd H = MatrixXd::Zero(SUNFLOWER_NR, STATE_SIZE);
    for (int i = 0; i < SUNFLOWER_NR; i++) {
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

TrilatMeasurement TrilatEKF::toTrilatMeasurement(Measurement m0, Measurement m1, Measurement m2){
    TrilatMeasurement tm;
    tm.timestamp_ = m0.timestamp_;
    tm.sensorLocs_ = MatrixXd(SUNFLOWER_NR, 2);
    tm.sensorLocs_ << m0.sensorLoc_.transpose(), m1.sensorLoc_.transpose(), m2.sensorLoc_.transpose();
    tm.distances_ = VectorXd(SUNFLOWER_NR);
    tm.distances_ << m0.distance_, m1.distance_, m2.distance_;
    //std::cout << "locs: " << tm.sensorLocs_ << std::endl;
    //std::cout << "dist: " << tm.distances_ << std::endl;
    return tm;
}

vector<TrilatMeasurement> TrilatEKF::getCombinations(vector<Measurement> mVec){
    // TODO make this more versatile
    vector<TrilatMeasurement> tmVec(8);
    tmVec[0] = toTrilatMeasurement(mVec[0], mVec[2], mVec[4]);
    tmVec[1] = toTrilatMeasurement(mVec[0], mVec[2], mVec[5]);
    tmVec[2] = toTrilatMeasurement(mVec[0], mVec[3], mVec[4]);
    tmVec[3] = toTrilatMeasurement(mVec[0], mVec[3], mVec[5]);
    tmVec[4] = toTrilatMeasurement(mVec[1], mVec[3], mVec[5]);
    tmVec[5] = toTrilatMeasurement(mVec[1], mVec[3], mVec[4]);
    tmVec[6] = toTrilatMeasurement(mVec[1], mVec[2], mVec[5]);
    tmVec[7] = toTrilatMeasurement(mVec[1], mVec[2], mVec[4]);
    
    return tmVec;
}

