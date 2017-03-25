#!/bin/bash

g++ -std=gnu++11 -I /home/sandipan/Point_process/dlib-19.4/ -o point_mle mle_point_process.cpp optimization_function.h optimization_function.cpp
