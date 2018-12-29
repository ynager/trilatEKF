/** @file kalman.cpp
 *  @brief Kalman filter function implementations
 *  @author Yannik Nager
 */

#include <math.h>
#include <iostream>
#include "kalman.h"

Kalman::Kalman() {}

Kalman::~Kalman() {}

void Kalman::initialize(VectorXd &x, MatrixXd &P, MatrixXd &F,
                        MatrixXd &H, MatrixXd &R, MatrixXd &Q) {
    x_ = x;
    P_ = P;
    F_ = F;
    H_ = H;
    R_ = R;
    Q_ = Q;
}

void Kalman::predict() {
    x_ = F_ * x_;   // project state
    P_ = F_ * P_ * F_.transpose() + Q_; // project error covariance
    //std::cout << P_ << std::endl;
    //std::cout << x_ << std::endl;
}

void Kalman::update(const VectorXd &z) {
    //std::cout << "z: " << z << std::endl;
    //std::cout << "x: " << x_ << std::endl;
    VectorXd y = z - H_ * x_;       // D - H*X_k (innovation)
    
    MatrixXd H_t = H_.transpose();  // transpose measurement matrix
    MatrixXd PH_t = P_ * H_t;
    MatrixXd S = H_ * PH_t + R_;    // innovation covariance
    MatrixXd S_i = S.inverse();
    MatrixXd K = PH_t * S_i;        // kalman gain
    
    update(y, K);                   // perform updates
}

uint16_t Kalman::updateMahalanobis(const MatrixXd zVec) {
    
    double dist_min = std::numeric_limits<double>::max();
    uint16_t idx_best = 0;
    VectorXd y_best;
    MatrixXd PH_t_best;
    MatrixXd S_i_best;
    
    for(int i = 0; i < zVec.cols(); ++i) {
        
        VectorXd z = zVec.col(i);       // extract single measurement
        VectorXd y = z - H_ * x_;       // D - H*X_k (innovation)
        
        MatrixXd H_t = H_.transpose();  // transpose measurement matrix
        MatrixXd PH_t = P_ * H_t;
        MatrixXd S = H_ * PH_t + R_;    // innovation covariance
        MatrixXd S_i = S.inverse();
        
        // calculate Mahalanobis distance
        double dist = z.transpose() * S_i * z;
        if(dist < dist_min) {
            idx_best = i;
            y_best = y;
            PH_t_best = PH_t;
            S_i_best = S_i;
            dist_min = dist;
        }
        
    }
    MatrixXd K = PH_t_best * S_i_best;  // kalman gain
    update(y_best, K);                  // perform updates
    
    return idx_best;
}

void Kalman::update(const VectorXd &y, const MatrixXd &K) {
    x_ = x_ + (K * y);             // update x_
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;             // update p_
}




