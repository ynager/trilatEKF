#include <math.h>
#include "kalman.h"

#define SUNFLOWER_NR 3


Kalman::Kalman(VectorXd &x, MatrixXd &P, MatrixXd &F,
               MatrixXd &H, MatrixXd &R, MatrixXd &Q) :
x_(x),
P_(P),
F_(F),
Q_(Q),
H_(H),
R_(R) {}


Kalman::~Kalman() {}

void Kalman::predict() {
    x_ = F_ * x_;   // project state
    MatrixXd F_t = F_.transpose();
    P_ = F_ * P_ * F_t + Q_; // project error covariance
}

void Kalman::update(const VectorXd &z) {
    VectorXd y = z - x_; //
    
    // build measurement matrix H using trilateration of measurements
    
    
    // get new measurement matrix from jacobian
    H_ = trilaterationJacobian(measurement);
    
    // standart equations
    MatrixXd H_t = H_.transpose(); // transpose measurement matrix
    MatrixXd PH_t = P_ * H_t;
    MatrixXd S = H_ * PH_t + R_; // innovation covariance
    MatrixXd S_i = S.inverse();
    MatrixXd K = PH_t * S_i; // kalman gain
    
    // calculate new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

MatrixXd Kalman::trilaterationJacobian(const VectorXd &z, const MatrixXd &beaconLoc){
    MatrixXd H = MatrixXd::Zero(SUNFLOWER_NR, 4);
    
    for (int i = 0; i < SUNFLOWER_NR; i++) {
        double den = sqrt(pow((x_(0)-beaconLoc(i,0)),2) + pow(x_(1)-beaconLoc(i,1),2));
        double dh_dx = (x_(0)-beaconLoc(i,1)) / den;
        double dh_dy = (x_(1)-beaconLoc(i,2)) / den;
        H(i, 0) = dh_dx;
        H(i, 1) = dh_dy;
    }
    return H;
}
