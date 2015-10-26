/*
This model is the Inverted Pendulum problem found in the paper
"lease-squared policy iterations"
*/
#pragma once
#include <string>
#include <cmath>
#include <iostream>
#include <random>
#include "common_RL.h"
#include "BaseModel.h"

namespace ReinforcementLearning {  
class Model_QuadrupoleBase: public BaseModel{
public:
    Model_QuadrupoleBase(){}
    virtual ~Model_QuadrupoleBase(){}
    virtual void run(int action) = 0;
    virtual void createInitialState() = 0;
    virtual double getRewards() = 0;
    virtual bool terminate() = 0;
protected:
    double dt;
    int trajOutputInterval;
    int timeCounter,fileCounter;
    double reward;
    std::random_device rand_generator;
    std::mt19937 rand_gen_mt19937;
    std::shared_ptr<std::normal_distribution<double>> rand_normal;
};
}
