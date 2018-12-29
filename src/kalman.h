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
     * @param x_ Initial state
     * @param P_ Initial state covariance
     * @param F_ Transistion matrix
     * @param H_ Measurement matrix
     * @param R_ Measurement covariance matrix
     * @param Q_ Process covariance matrix
     */
    void initialize(VectorXd &x_, MatrixXd &P_, MatrixXd &F_, MatrixXd &H_, MatrixXd &R_, MatrixXd &Q_);
    
    /**
     * @brief Predict state and state covariance
     */
    void predict();
    
    /**
     * @brief Update state and state covariance using extended kalman filter
     * @param z_ measurement at timestep k+1
     */
    void update(const VectorXd &z_);
    
    /**
     * @brief Update state and state covariance using measurement with lowest Mahalanobis distance
     * @param zVec_ vector of measurements at timestep k+1
     * @return idx of measurement with lowest Mahalanobis distance
     */
    uint16_t updateMahalanobis(const MatrixXd zVec);
    
private:
    
    /**
     * @brief Update state and state covariance using extended kalman filter
     * @param y
     * @param K kalman gain
     */
    void update(const VectorXd &y, const MatrixXd &K);
};

#endif //KALMAN_H_
