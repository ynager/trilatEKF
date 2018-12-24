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
#define STATE_SIZE 4
#define VAR_Z 0.1 // sunflower meas variance
#define MEAS_NOISE_X 5.05
#define MEAS_NOISE_Y 5.05

using namespace Eigen;

/** Single Sunflower measurement */
struct Measurement {
    long long timestamp_;
    Vector2d sensorLoc_;
    double distance_;
};

/** Trilateration Measurement struct */
struct TrilatMeasurement {
    long long timestamp_;
    MatrixXd sensorLocs_;
    VectorXd distances_;
};

class TrilatEKF{
public:
    
    /** Kalman Filter */
    Kalman ekf_;
    
    /**
     * Constructor
     */
    TrilatEKF(VectorXd xInit, MatrixXd sensorLoc);
    
    /**
     * Destructor
     */
    ~TrilatEKF();
    
    /**
     * Process measurement of type TrilatMeasurement
     * @param Trilateration measurement package
     */
    void processMeasurement(const TrilatMeasurement &trilatMeasurement);
    
    /**
     * calculate jacobian of measurement matrix
     * @param measurement of type TrilatMeasurement
     */
    MatrixXd getJacobian(const MatrixXd &sensorLoc, const VectorXd &x);
    
    /**
     * Match single measurements into trilat measurement
     */
    TrilatMeasurement matchMeasurements(std::vector<Measurement> mvec);
    
private:
    
    Eigen::MatrixXd R_;
    Eigen::MatrixXd H_;
    long long timestamp_prev_;
    
};

#endif
