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

#define STATE_SIZE 4

using namespace std;
using namespace boost;
using Eigen::MatrixXd;
using Eigen::VectorXd;

int main(int argc, char* argv[]) {
    cout << "Running trilatEKF...\n";
    
    // set EKF params
    TrilatParams ekfparams = {0.1, 10.0, 10.0};
    
    // initialize objects
    VectorXd x_init_0(STATE_SIZE);
    VectorXd x_init_1(STATE_SIZE);
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
    //string data("../data/dataset_simple2.csv"); // load dataset
    
    // Input file handling
    fstream in(data.c_str());
    if (!in.is_open()) return 1;
    typedef tokenizer< escaped_list_separator<char> > Tokenizer;
    vector<string> vec;
    string line;
    long cnt = -1;
    std::vector<Measurement> mVec;
    
    // Loop through measurements
    while (getline(in,line) )
    {
        Tokenizer tok(line);
        vec.assign(tok.begin(),tok.end());
        
        if (cnt == -1) {            // set initial position estimates
            x_init_1 << stod(vec[0]), stod(vec[1]), 0.0, 0.0; //, 0.0, 0.0;
            x_init_0 << stod(vec[2]), stod(vec[3]), 0.0, 0.0; //, 0.0, 0.0;
            tekf_0 = new TrilatEKF(x_init_0, sensorloc, ekfparams);
            tekf_1 = new TrilatEKF(x_init_1, sensorloc, ekfparams);
        }
        else {
            // assign to measurement struct
            Measurement m;
            m.timestamp = stoi(vec[0]);
            m.sensorloc << stod(vec[1]), stod(vec[2]);
            m.distance = stod(vec[3]);
            mVec.push_back(m);
            
            if (cnt % 6 == 5) { // after 6 measurements generate combinations and run EKF
                // generate all measurement combinations
                std::vector<TrilatMeasurement> tmVec = tekf_0->getCombinations(mVec);
                mVec.clear();
                tekf_0->processMeasurements(&tmVec);
                tekf_1->processMeasurements(&tmVec);
                
                // console outputs
                std::cout << "Timestamp: " << m.timestamp << std::endl;
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
            }
        }
        cnt += 1;
    }
    outfile.close(); // close output file
}
