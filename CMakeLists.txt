
cmake_minimum_required(VERSION 3.0.0)
project(trilatEKF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(include ${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})


set(SOURCE_FILES src/kalman.cpp src/main.cpp src/trilatEKF.cpp)
add_executable(main ${SOURCE_FILES})
