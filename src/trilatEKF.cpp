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
    
    TrilatMeasurement m0 = (*measurements)[0];
    double dt = (m0.timestamp - timestamp_prev_) / 1000.0;  // time delta [seconds]
    if(dt == 0) return;
    double dt2 = dt * dt;
    double dt3 = dt2 * dt;
    double dt4 = dt3 * dt;
    double dt5 = dt4 * dt;
    timestamp_prev_ = m0.timestamp;
    
    // PREDICTION
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
            ekf_.Q_ <<  p_.var_x, 0, 0, p_.var_y;
    }
    
    ekf_.predict();     // perform prediction step
    
    // UPDATE
    MatrixXd H = getJacobian(m0.sensorlocs, ekf_.x_);   // get current jacobian
    ekf_.H_ = H;    // TODO set H row to zero if sensor out of range
    
    // put distances into Matrix
    MatrixXd distanceMatrix(SUNFLOWER_NR, (*measurements).size());
    for (int i = 0; i < (*measurements).size(); i++)
        distanceMatrix.col(i) = (*measurements)[i].distances;
    
    // calculate expected measurement
    VectorXd h(3);
    h <<    sqrt( pow(ekf_.x_(0) - m0.sensorlocs(0,0), 2) + pow(ekf_.x_(1) - m0.sensorlocs(0, 1), 2)),
            sqrt( pow(ekf_.x_(0) - m0.sensorlocs(1,0), 2) + pow(ekf_.x_(1) - m0.sensorlocs(1, 1), 2)),
            sqrt( pow(ekf_.x_(0) - m0.sensorlocs(2,0), 2) + pow(ekf_.x_(1) - m0.sensorlocs(2, 1), 2));
    
    // perform update step
    uint16_t idx = ekf_.updateMahalanobis(distanceMatrix, h);
    
    // remove "used" measurements
    // measurements for other object known since only two readings per sensor
    TrilatMeasurement a = (*measurements)[(idx + 4) % 8];
    (*measurements).clear();
    (*measurements).push_back(a);
}

MatrixXd TrilatEKF::getJacobian(const MatrixXd &sensorlocs, const VectorXd &x){
    MatrixXd H = MatrixXd::Zero(SUNFLOWER_NR, statesize_);
    for (int i = 0; i < SUNFLOWER_NR; ++i) {
        double den = sqrt(pow(x(0)-sensorlocs(i,0),2) + pow(x(1)-sensorlocs(i,1),2));
        if (den < 0.001) {
            std::cout << "Error: Division by zero in getJacobian()\n";
            return H;
        }
        double dh_dx = (x(0)-sensorlocs(i,0)) / den;
        double dh_dy = (x(1)-sensorlocs(i,1)) / den;
        H(i, 0) = dh_dx;
        H(i, 1) = dh_dy;
    }
    return H;
}

vector<TrilatMeasurement> TrilatEKF::getCombinations(const vector<Measurement> &mvec){
    assert(mvec.size() == 6); // TODO make this more versatile
    vector<TrilatMeasurement> tmvec(8);
    tmvec[0] = toTrilatMeasurement(mvec[0], mvec[2], mvec[4]);
    tmvec[1] = toTrilatMeasurement(mvec[0], mvec[2], mvec[5]);
    tmvec[2] = toTrilatMeasurement(mvec[0], mvec[3], mvec[4]);
    tmvec[3] = toTrilatMeasurement(mvec[0], mvec[3], mvec[5]);
    tmvec[4] = toTrilatMeasurement(mvec[1], mvec[3], mvec[5]);
    tmvec[5] = toTrilatMeasurement(mvec[1], mvec[3], mvec[4]);
    tmvec[6] = toTrilatMeasurement(mvec[1], mvec[2], mvec[5]);
    tmvec[7] = toTrilatMeasurement(mvec[1], mvec[2], mvec[4]);
    return tmvec;
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
