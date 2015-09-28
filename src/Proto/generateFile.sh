#!/bin/bash

protoc ReinforcementLearning.proto --cpp_out=.
mv ReinforcementLearning.pb.h ../../include
mv ReinforcementLearning.pb.cc ../
