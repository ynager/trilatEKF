/** @file main.cpp
 *  @brief Main file running trilatEKF on an example dataset
 *  @author Yannik Nager
 */
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
// External
#include "Eigen/Dense"
#include "boost/tokenizer.hpp"
// Custom
#include "kalman.h"
#include "helpers.h"
#include "trilatEKF.h"

// set state size (2: P, 4: PV, 6: PVA)
#define STATE_SIZE 4

using namespace std;
using namespace boost;
using Eigen::MatrixXd;
using Eigen::VectorXd;

int main(int argc, char* argv[]) {
    cout << "Running trilatEKF...\n";
    
    // set EKF params (variance z, variance x, variance y)
    TrilatParams ekfparams = {0.01, 0.1, 0.1};
    
    // initialize objects
    VectorXd x_init_0 = VectorXd::Zero(STATE_SIZE);
    VectorXd x_init_1 = VectorXd::Zero(STATE_SIZE);
    MatrixXd sensorloc = MatrixXd(3,2); // initial sensor location
    sensorloc <<    0.0, 0.0,
                    5.0, 0.0,
                    0.0, 5.0;
    
    // Kalman filter for two objects
    TrilatEKF *tekf_0;
    TrilatEKF *tekf_1;
    
    // Output file handling
    std::ofstream outfile;
    outfile.open ("../data/output.csv");
    string data("../data/dataset_3.csv"); // load dataset
    
    // Input file handling
    fstream in(data.c_str());
    if (!in.is_open()) {
        cout << "Error: input file not found\n";
        return 1;
    }
    typedef tokenizer< escaped_list_separator<char> > Tokenizer;
    vector<string> vec;
    string line;
    
    // variable definitions
    Measurement m;
    std::vector<Measurement> mvec;
    long long timestamp_prev = 0;
    bool init = true;
    
    // Loop through measurements
    while (getline(in,line) )
    {
        Tokenizer tok(line);
        vec.assign(tok.begin(),tok.end());
        
        if (init) {    // set initial position estimates
            x_init_1(0) = stod(vec[0]);
            x_init_1(1) = stod(vec[1]);
            x_init_0(0) = stod(vec[2]);
            x_init_0(1) = stod(vec[3]);
            tekf_0 = new TrilatEKF(x_init_0, sensorloc, ekfparams);
            tekf_1 = new TrilatEKF(x_init_1, sensorloc, ekfparams);
            init = false;
        }
        else {
            long long timestamp = stoi(vec[0]);
            if (timestamp != timestamp_prev) { // check if timestamp has changed
                // generate all measurement combinations
                std::vector<TrilatMeasurement> tmvec = TrilatEKF::getCombinations(mvec);
                mvec.clear();
                
                // process measurement
                tekf_0->processMeasurements(&tmvec);
                tekf_1->processMeasurements(&tmvec);
                
                // console outputs
                std::cout << "Timestamp: " << timestamp << std::endl;
                std::cout << "Object 0 :" << tekf_0->ekf_.x_.transpose() << "\n";
                std::cout << "Object 1 :" << tekf_1->ekf_.x_.transpose() << "\n\n";
                
                // write to file
                outfile << m.timestamp << ","
                << 0 << ","
                << tekf_0->ekf_.x_(0) << ","
                << tekf_0->ekf_.x_(1) << "\n";
                
                outfile << m.timestamp << ","
                << 1 << ","
                << tekf_1->ekf_.x_(0) << ","
                << tekf_1->ekf_.x_(1) << "\n";
                
                timestamp_prev = timestamp;
            }
            // assign to measurement struct
            m.timestamp = stoi(vec[0]);
            m.sensorloc << stod(vec[1]), stod(vec[2]);
            m.distance = stod(vec[3]);
            mvec.push_back(m);
        }
    }
    outfile.close(); // close output file
}
