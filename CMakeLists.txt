
cmake_minimum_required(VERSION 3.0.0)
project(sunflower_EKF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(EIGEN_INC_DIR /usr/local/Cellar/eigen/3.3.7/include/eigen3)
include_directories(${EIGEN_INC_DIR})
include_directories(include)


set(SOURCE_FILES src/kalman.cpp src/main.cpp src/trilatEKF.cpp)
add_executable(main ${SOURCE_FILES})
