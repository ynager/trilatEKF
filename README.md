# TrilatEKF Project

This project aims to track moving objects through trilateration of distance measurements using an extended Kalman filter.

## Build instructions

This project requires [Eigen](http://eigen.tuxfamily.org) and [Boost](https://www.boost.org). To clone and build the project, type

```bash
git clone https://github.com/ynager/trilatEKF
cd trilatEKF
mkdir build
cd build
cmake ..
make
```

To run the executable type 
```bash
./trilatEKF
```

The two object state estimates are written to 'data/output.csv'


