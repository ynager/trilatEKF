/** @file trilatEKF.h
 *  @brief trilatEKF structs and class definition
 *  @author Yannik Nager
 */

#ifndef TRILAT_EKF_H_
#define TRILAT_EKF_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "Eigen/Dense"
#include "kalman.h"
#include "helpers.h"

#define SUNFLOWER_NR 3
#define VAR_Z 0.1 // sunflower meas variance
#define MEAS_NOISE_X 10.0
#define MEAS_NOISE_Y 10.0

using namespace Eigen;

/** @brief Single measurement */
struct Measurement {
    long long timestamp;
    Vector2d sensorloc;
    double distance;
};

/** Trilateration measurement struct */
struct TrilatMeasurement {
    long long timestamp;
    MatrixXd sensorlocs;
    VectorXd distances;
};

/** Trilateration EKF params */
struct TrilatParams {
    float var_z;
    float sigma_x;
    float sigma_y;
};

/**
 * TrilatEKF class
 */
class TrilatEKF{
public:
    /** Kalman Filter */
    Kalman ekf_;
    
    /**
     * @brief Constructor
     */
    TrilatEKF(VectorXd xInit, MatrixXd sensorLoc, TrilatParams p);
    
    /**
     * @brief Destructor
     */
    ~TrilatEKF();
    
    /**
     * @brief Process measurement of type TrilatMeasurement
     * @param Trilateration measurement package
     * @return void
     */
    void processMeasurements(std::vector<TrilatMeasurement>* trilatMeasurements);
    
    /**
     * @brief calculate jacobian of measurement matrix
     * @param sensorLoc current sensor location
     * @param x current state estimate
     * @return jacobian matrix
     */
    MatrixXd getJacobian(const MatrixXd &sensorLoc, const VectorXd &x);
    
    /**
     * @brief calculate jacobian of measurement matrix
     * @param m TrilatMeasurement
     * @return jacobian matrix
     */
    MatrixXd getJacobian(const TrilatMeasurement &m);
    
    /**
     * @brief Match single measurements into trilat measurement
     * @param m0 first measurement
     * @param m1 second measurement
     * @param m2 third measurement
     * @return TrilatMeasurement
     */
    TrilatMeasurement toTrilatMeasurement(Measurement m0, Measurement m1, Measurement m2);
    
    /**
     * @brief get measurement combinations
     * @param vVec vector containing 6 measurements (2 measurements per sensor)
     * @note measurements must be sorted by sensor
     * @return vector of TrilatMeasurements
     */
    std::vector<TrilatMeasurement> getCombinations(std::vector<Measurement> mvec);
    
private:
    
    Eigen::MatrixXd R_;
    Eigen::MatrixXd H_;
    Eigen::MatrixXd sensorloc_;
    TrilatParams p_;
    long long timestamp_prev_;
    int statesize_;
};

#endif
