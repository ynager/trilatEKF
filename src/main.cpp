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

using namespace std;
using namespace boost;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

int main(int argc, char* argv[]) {
    cout << "Sunflower EKF\n";
    
    // initialize objects
    VectorXd xInit_a(STATE_SIZE);
    VectorXd xInit_b(STATE_SIZE);
    MatrixXd sensorLoc = MatrixXd(3,2); // initial sensor location
    sensorLoc << 0.0, 0.0,
    5.0, 0.0,
    0.0, 5.0;
    
    // Kalman Filter for objecs
    TrilatEKF *tEKF_0;
    TrilatEKF *tEKF_1;
    
    // Output file handling
    std::ofstream myfile;
    myfile.open ("output.csv");
    
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
    
    while (getline(in,line) ) //&& cnt < 20)
    {
        Tokenizer tok(line);
        vec.assign(tok.begin(),tok.end());
        
        if (cnt == -1) {            // set initial position estimates
            xInit_a << stod(vec[0]), stod(vec[1]), 0.0, 0.0; //, 0.0, 0.0;
            xInit_b << stod(vec[2]), stod(vec[3]), 0.0, 0.0; //, 0.0, 0.0;
            tEKF_0 = new TrilatEKF(xInit_a, sensorLoc);
            tEKF_1 = new TrilatEKF(xInit_b, sensorLoc);
        }
        else {
            // assign to measurement struct
            Measurement m;
            m.timestamp_ = stoi(vec[0]);
            m.sensorLoc_ << stod(vec[1]), stod(vec[2]);
            m.distance_ = stod(vec[3]);
            mVec.push_back(m);
            
            if (cnt % 6 == 5) { // after 6 measurements match & run EKF
                // generate all measurement combinations
                std::vector<TrilatMeasurement> tmVec = tEKF_0->getCombinations(mVec);
                mVec.clear();
                tEKF_0->processMeasurements(tmVec);
                tEKF_1->processMeasurements(tmVec);
                
                // write to file
                myfile << m.timestamp_ << ","
                << 0 << ","
                << tEKF_0->ekf_.x_(0) << ","
                << tEKF_0->ekf_.x_(1) << "\n";
                
                myfile << m.timestamp_ << ","
                << 1 << ","
                << tEKF_1->ekf_.x_(0) << ","
                << tEKF_1->ekf_.x_(1) << "\n";
            }
        }
        cnt += 1;
    }
    myfile.close(); // close output file
}
