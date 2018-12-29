/** @file trilatEKF.cpp
 *  @brief trilatEKF method implementations
 *  @author Yannik Nager
 */

#include "trilatEKF.h"
using namespace std;

TrilatEKF::TrilatEKF(const VectorXd &x_init, const MatrixXd &sensorloc, const TrilatParams &p) {
    
    // TODO improve initialization
    this->statesize_ = x_init.size();
    this->p_ = p;
    this->timestamp_prev_ = 0;
    this->sensorloc_ = sensorloc;
    
    MatrixXd F_init = MatrixXd::Identity(statesize_, statesize_);                   // set state transition matrix F
    MatrixXd Q_init = MatrixXd::Identity(statesize_, statesize_);                   // set process noise covariance
    
    MatrixXd H_init = getJacobian(sensorloc, x_init);                               // calculate initial measurement matrix H
    MatrixXd R_init = MatrixXd::Identity(SUNFLOWER_NR, SUNFLOWER_NR) * p_.var_z;    // set measurement covariance matrix R
    MatrixXd P_init = MatrixXd::Identity(statesize_, statesize_)*10;                // only position of initial state is known
    P_init(0,0) = 0.01; // TODO put outside TrilatEKF
    P_init(1,1) = 0.01;
    
    ekf_.initialize(x_init, P_init, F_init, H_init, R_init, Q_init);                // initialize Kalman Filter with x_init
}

TrilatEKF::~TrilatEKF() {}

void TrilatEKF::processMeasurements(std::vector<TrilatMeasurement>* measurements){
    
    // PREDICTION
    double dt = ((*measurements)[0].timestamp - timestamp_prev_) / 1000.0;  // get time delta dt [seconds]
    
    if(dt == 0) {
        return;
    }
    
    double dt2 = dt * dt;
    double dt3 = dt2 * dt;
    double dt4 = dt3 * dt;
    double dt5 = dt4 * dt;
    timestamp_prev_ = (*measurements)[0].timestamp;
    
    // update F and Q matrices
    switch(statesize_) {
        case 6:
            ekf_.F_ <<  1, 0, dt, 0, dt2/2, 0,
            0, 1, 0, dt, 0, dt2/2,
            0, 0, 1, 0, dt, 0,
            0, 0, 0, 1, 0, dt,
            0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 1;
            ekf_.Q_ <<  dt5/20*p_.var_x, 0, dt4/8*p_.var_x, 0, dt3/6*p_.var_x, 0,
            0, dt5/20*p_.var_y, 0, dt4/8*p_.var_y, 0, dt3/6*p_.var_y,
            dt4/8*p_.var_x, 0, dt3/3*p_.var_x, 0, dt2/2*p_.var_x, 0,
            0, dt4/8*p_.var_y, 0, dt3/3*p_.var_y, 0, dt2/2*p_.var_x,
            dt3/3*p_.var_x, 0, dt2/2*p_.var_x, 0, dt*p_.var_x, 0 ,
            0, dt3/3*p_.var_y, 0, dt2/2*p_.var_y, 0, dt*p_.var_y;
            break;
            
        case 4:
            ekf_.F_(0, 2) = dt;
            ekf_.F_(1, 3) = dt;
            ekf_.Q_ <<  dt4/4*p_.var_x, 0, dt3/2*p_.var_x, 0,
            0, dt4/4*p_.var_y, 0, dt3/2*p_.var_y,
            dt3/2*p_.var_x, 0, dt2*p_.var_x, 0,
            0, dt3/2*p_.var_y, 0, dt2*p_.var_y;
            break;
            
        case 2:
            ekf_.Q_ << p_.var_x, 0, 0, p_.var_y;
    }
    
    ekf_.predict();     // perform prediction step
    
    // UPDATE
    MatrixXd H = getJacobian((*measurements)[0].sensorlocs, ekf_.x_);   // get current jacobian
    ekf_.H_ = H;
    // TODO set H row to zero if sensor out of range
    
    // put distances into Matrix
    MatrixXd distanceMatrix(SUNFLOWER_NR, (*measurements).size());
    for (int i = 0; i < (*measurements).size(); i++)
        distanceMatrix.col(i) = (*measurements)[i].distances;
    
    // calculate expected measurement
    VectorXd h(3);
    h <<    sqrt( pow(ekf_.x_(0) - sensorloc_(0,0), 2) + pow(ekf_.x_(1) - sensorloc_(0, 1), 2)),
            sqrt( pow(ekf_.x_(0) - sensorloc_(1,0), 2) + pow(ekf_.x_(1) - sensorloc_(1, 1), 2)),
            sqrt( pow(ekf_.x_(0) - sensorloc_(2,0), 2) + pow(ekf_.x_(1) - sensorloc_(2, 1), 2));
    
    // perform update step
    uint16_t idx = ekf_.updateMahalanobis(distanceMatrix, h);
    
    // remove "used" measurements
    TrilatMeasurement a = (*measurements)[(idx + 4) % 8]; // measurements for other object known since only two readings per sensor
    (*measurements).clear();
    (*measurements).push_back(a);
}

MatrixXd TrilatEKF::getJacobian(const MatrixXd &sensorloc, const VectorXd &x){
    MatrixXd H = MatrixXd::Zero(SUNFLOWER_NR, statesize_);
    for (int i = 0; i < SUNFLOWER_NR; ++i) {
        double den = sqrt(pow(x(0)-sensorloc(i,0),2) + pow(x(1)-sensorloc(i,1),2));
        if (den < 0.001) {
            std::cout << "Error: Division by zero in getJacobian()\n";
            return H;
        }
        double dh_dx = (x(0)-sensorloc(i,0)) / den;
        double dh_dy = (x(1)-sensorloc(i,1)) / den;
        H(i, 0) = dh_dx;
        H(i, 1) = dh_dy;
    }
    return H;
}

vector<TrilatMeasurement> TrilatEKF::getCombinations(const vector<Measurement> &mVec){
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

TrilatMeasurement TrilatEKF::toTrilatMeasurement(const Measurement &m0, const Measurement &m1, const Measurement &m2){
    TrilatMeasurement tm;
    tm.timestamp = m0.timestamp;
    tm.sensorlocs = MatrixXd(SUNFLOWER_NR, 2);
    tm.sensorlocs << m0.sensorloc.transpose(), m1.sensorloc.transpose(), m2.sensorloc.transpose();
    tm.distances = VectorXd(SUNFLOWER_NR);
    tm.distances << m0.distance, m1.distance, m2.distance;
    return tm;
}

