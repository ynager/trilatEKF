#ifndef KALMAN_H_
#define KALMAN_H_

#include "Eigen/Dense"

using namespace Eigen;

class Kalman {
public:
    
    VectorXd x_;  // state
    MatrixXd P_;  // covariance matrix
    MatrixXd F_;  // state transistion matrix
    MatrixXd Q_;  // process covariance matrix
    MatrixXd H_;  // measurement matrix
    MatrixXd R_;  // measurement covariance matrix
    
    /** Constructor - Initialize Kalman Filter
     @param x_ Initial state
     @param P_ Initial state covariance
     @param F_ Transistion matrix
     @param H_ Measurement matrix
     @param R_ Measurement covariance matrix
     @param Q_ Process covariance matrix
     */
    Kalman(VectorXd &x_, MatrixXd &P_, MatrixXd &F_, MatrixXd &H_, MatrixXd &R_, MatrixXd &Q_);
    
    /* Destructor */
    ~Kalman();
    
    /**
     Predict state and state covariance
     */
    void predict();
    
    /**
     Update state and state covariance using extended kalman filter
     @param z_ measurement at timestep k+1
     */
    void update(const VectorXd &z_);
    
private:
    
    MatrixXd trilaterationJacobian(const VectorXd &z, const MatrixXd &beaconLoc);
    
};

#endif //KALMAN_H_
