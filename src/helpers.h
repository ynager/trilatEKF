#ifndef HELPERS_H_
#define HELPERS_H_

#include "Eigen/Dense"

/** Sunflower Measurement struct */
struct FlowerMeasurement {
    long long timestamp_;
    Eigen::MatrixXd sensorLocs_;
    Eigen::VectorXd distances_;
};

#endif
