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

/** Sunflower Measurement struct */
struct TrilatMeasurement {
    long long timestamp_;
    MatrixXd sensorLoc_;
    VectorXd distance_;
};

class TrilatEKF{
public:
    
    Kalman ekf_;
    
    /** Constructor */
    TrilatEKF(VectorXd xInit, MatrixXd sensorLoc);
    
    /** Destructor */
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
    
private:
    
    Eigen::MatrixXd R_;
    Eigen::MatrixXd H_;
    long long timestamp_prev_ = 0;
    
};

#endif
