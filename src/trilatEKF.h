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
#define MEAS_NOISE_X 10.0
#define MEAS_NOISE_Y 10.0

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
    void processMeasurements(const std::vector<TrilatMeasurement> &trilatMeasurements);
    
    /**
     * calculate jacobian of measurement matrix
     * @param measurement of type TrilatMeasurement
     */
    MatrixXd getJacobian(const MatrixXd &sensorLoc, const VectorXd &x);
    
    MatrixXd getJacobian(const TrilatMeasurement &m);
    
    /**
     * Match single measurements into trilat measurement
     */
    TrilatMeasurement toTrilatMeasurement(Measurement m0, Measurement m1, Measurement m2);
    
    std::vector<TrilatMeasurement> getCombinations(std::vector<Measurement> mVec);
    
private:
    
    Eigen::MatrixXd R_;
    Eigen::MatrixXd H_;
    Eigen::MatrixXd sensorLoc_;
    long long timestamp_prev_;
};

#endif
