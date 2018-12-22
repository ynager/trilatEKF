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
#define DT 0.25 // TODO  check parameters
#define ACCEL_NOISE_MAG 0.000001/DT
#define VAR_Z 0.1 // sunflower meas variance
#define MEAS_NOISE_X 3.0
#define MEAS_NOISE_Y 3.0

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
    TrilatMeasurement matchMeasurements(std::vector<Measurement> mvec){
        TrilatMeasurement tm;
        tm.timestamp_ = mvec[0].timestamp_;
        tm.sensorLocs_ = MatrixXd(2,3); // TODO magic var
        tm.distances_ = VectorXd(3);
        tm.distances_ << mvec[0].distance_, mvec[0].distance_, mvec[0].distance_;
        return tm;
    }
    
private:
    
    Eigen::MatrixXd R_;
    Eigen::MatrixXd H_;
    long long timestamp_prev_ = 0;
    
};

#endif
