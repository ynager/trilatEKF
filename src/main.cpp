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
    Vector4d xInit_a;
    Vector4d xInit_b;
    MatrixXd sensorLoc = MatrixXd(3,2); // initial sensor location
    sensorLoc << 0.0, 0.0,
                 5.0, 0.0,
                 0.0, 5.0;
    
    // Trilateration Kalman Filter for object A
    TrilatEKF *tEKF;
    
    //string data("../data/dataset_3.csv"); // load dataset
    string data("../data/dataset_simple2.csv"); // load dataset

    fstream in(data.c_str());
    if (!in.is_open()) return 1;
    
    typedef tokenizer< escaped_list_separator<char> > Tokenizer;
    vector<string> vec;
    string line;
    
    
    long cnt = -1;
    // process dataset line by line
    std::vector<Measurement> mvec(3);
    
    // output file
    std::ofstream myfile;
    myfile.open ("output.csv");
    
    while (getline(in,line))
    {
        Tokenizer tok(line);
        vec.assign(tok.begin(),tok.end());
        //copy(vec.begin(), vec.end(), ostream_iterator<string>(cout, " "));
        
        // set initial position estimates
        if (cnt == -1) {
            xInit_a << stod(vec[0]), stod(vec[1]), 0.0, 0.0;
            xInit_b << stod(vec[2]), stod(vec[3]), 0.0, 0.0;
            tEKF = new TrilatEKF(xInit_a, sensorLoc);
        }
        else {
            // assign data to individual measurement structs
                Measurement m;
                m.timestamp_ = stoi(vec[0]);
                m.sensorLoc_ << stod(vec[1]), stod(vec[2]);
                m.distance_ = stod(vec[3]);
                mvec[(cnt%3)] = m; // fill in measurements
                
                //copy(vec.begin(), vec.end(), ostream_iterator<string>(cout, " "));
                
                if (cnt % 3 == 2) { // after 3 measurements run EKF
                    cout << "meas: " << mvec[0].distance_ << ", " << mvec[1].distance_ << ", " << mvec[2].distance_ << std::endl;
                    TrilatMeasurement tm = tEKF->matchMeasurements(mvec);
                    tEKF->processMeasurement(tm);
                    cout << "x: " << tEKF->ekf_.x_.transpose() << endl;
                    
                    myfile << tEKF->ekf_.x_(0) << ","
                           << tEKF->ekf_.x_(1) << ","
                           << tEKF->ekf_.x_(2) << ","
                           << tEKF->ekf_.x_(3) << "\n";

            }
        }
        cnt += 1;
    }
    myfile.close(); // close output file
}
