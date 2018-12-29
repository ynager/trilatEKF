/** @file kalman.h
 *  @brief Kalman filter class definition
 *  @author Yannik Nager
 */

#ifndef KALMAN_H_
#define KALMAN_H_

#include <vector>
#include "Eigen/Dense"

using namespace Eigen;

/**
 * Kalman filter class
 */
class Kalman {
public:
    
    VectorXd x_;    // state
    MatrixXd P_;    // covariance matrix
    MatrixXd F_;    // state transistion matrix
    MatrixXd Q_;    // process covariance matrix
    MatrixXd H_;    // measurement matrix
    MatrixXd R_;    // measurement covariance matrix
    // TODO: make access more restricted
    
    /**
     * @brief Constructor
     */
    Kalman();
    
    /**
     * @brief Destructor
     */
    ~Kalman();
    
    /**
     * @brief Initialize Kalman Filter
     * @param x Initial state estimate
     * @param P Initial state covariance
     * @param F Transistion matrix
     * @param H Measurement matrix
     * @param R Measurement covariance matrix
     * @param Q Process covariance matrix
     */
    void initialize(const VectorXd &x, const MatrixXd &P, const MatrixXd &F,
                    const MatrixXd &H, const MatrixXd &R, const MatrixXd &Q);
    
    /**
     * @brief Predict state and state covariance
     */
    void predict();
    
    /**
     * @brief Update state and state covariance using extended kalman filter
     * @param z measurement at timestep k+1
     */
    void update(const VectorXd &z);
    
    /**
     * @brief Update state and state covariance using measurement with lowest Mahalanobis distance
     * @param zvec vector of measurements at timestep k+1
     * @param h expected measurement
     * @return idx of measurement with lowest Mahalanobis distance
     */
    uint16_t updateMahalanobis(const MatrixXd &zvec, const VectorXd &h);
    
private:
    
    /**
     * @brief Update state and state covariance using extended kalman filter
     * @param y
     * @param K kalman gain
     */
    void update(const VectorXd &y, const MatrixXd &K);
};

#endif //KALMAN_H_
