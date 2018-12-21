#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
// External
#include "Eigen/Dense"
// Custom
#include "kalman.h"
#include "helpers.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

int main(int argc, char* argv[]) {
    cout << "Sunflower EKF\n";
    
    MatrixXd Q = MatrixXd::Zero(4, 4);
    cout << Q << endl;
    
}
