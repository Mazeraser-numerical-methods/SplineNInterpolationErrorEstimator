cmake_minimum_required(VERSION 4.3)
project(Spline_n_InterporationErrorEstimator)

set(CMAKE_CXX_SCAN_FOR_MODULES OFF CACHE BOOL "" FORCE)

set(CMAKE_CXX_STANDARD 20)

add_executable(Spline_n_InterporationErrorEstimator main.cpp
        point.h
        cubic_spline.cpp
        cubic_spline.h
        interpolation_error_estimator.cpp
        interpolation_error_estimator.h)
