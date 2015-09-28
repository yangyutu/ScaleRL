#pragma once
#include <string>
#include <cmath>
#include <iostream>
#include <random>
#include "Model_QuadrupoleBase.h"
#include "boost/multi_array.hpp"

namespace ReinforcementLearning {  
class Model_QuadrupoleLD: public Model_QuadrupoleBase{
public:
    Model_QuadrupoleLD(){}
    Model_QuadrupoleLD(std::string filetag0);
    virtual ~Model_QuadrupoleLD(){}
    virtual void run(int action);
    virtual void createInitialState();
    virtual double getRewards();
    virtual bool terminate();
protected:
    double xminop, yminop, xdelop, ydelop;
    int xbinNum, ybinNum;
    int trajOutputInterval;
    int timeCounter,fileCounter;
    bool trajOutputFlag;
    std::ofstream trajOs;
    std::string filetag;
    typedef boost::multi_array<double, 4> Array4D_type;
    typedef boost::multi_array<int, 2> Array2D_type;
    std::shared_ptr<Array4D_type> v, d;
    std::shared_ptr<Array2D_type> hasValue;
    void load_V_D(std::string filename[]);
    void outputTrajectory(std::ostream& os, int action);
};
}
