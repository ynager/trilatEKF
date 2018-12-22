#include <math.h>
#include <iostream>
#include "kalman.h"

#define SUNFLOWER_NR 3


Kalman::Kalman() {}

Kalman::~Kalman() {}

void Kalman::initialize(VectorXd &x, MatrixXd &P, MatrixXd &F,
                        MatrixXd &H, MatrixXd &R, MatrixXd &Q) {
    x_ = x;
    P_ = P;
    F_ = F;
    Q_ = Q;
    H_ = H;
    R_ = R;
}

void Kalman::predict() {
    //std::cout << F_ << std::endl;
    //std::cout << x_ << std::endl;
    x_ = F_ * x_;   // project state
    MatrixXd F_t = F_.transpose();
    P_ = F_ * P_ * F_t + Q_; // project error covariance
}

void Kalman::update(const VectorXd &z) {
    //std::cout << "z: " << z << std::endl;
    //std::cout << "x: " << x_ << std::endl;
    VectorXd y = z - H_ * x_; // D - H*X_k
    
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
