#pragma once
#include <string>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include "Model_QuadrupoleBase.h"
#include <math.h>
#include <armadillo>
#include <stdlib.h>

namespace ReinforcementLearning {  
class Model_QuadrupoleMC: public Model_QuadrupoleBase{
public:
    Model_QuadrupoleMC(){}
    Model_QuadrupoleMC(std::string filetag0, int R, int polygon);
    virtual ~Model_QuadrupoleMC(){ }
    virtual void run(int action);
    virtual void createInitialState();
    virtual double getRewards();
    virtual bool terminate();
protected:
    static const int np = 300;
    static const int np3 = 900;
    static const int IndexR = 60;
    const double pi = 3.1415925025939941;
    const int OverLapCheck = 0;
    std::vector<double> RCheck, Edge;
    arma::field<std::vector<int> > IndexMap;
    int R, n_rows, n_cols, polygon, ppp;
    double dx1, dx2;
    double Angle, EdgeLength, a;
    int DiscretizedR[np][2];
    int opt;
    double nstep;
    int trajOutputInterval;
    int timeCounter,fileCounter;
    bool trajOutputFlag;
    std::ofstream trajOs, opOs;
    std::string filetag;
    std::ofstream file;
    
    double r[np3], psi6, rg, lambda, rmin, dt;
    double DiffTrans, DiffRot;
    std::shared_ptr<std::uniform_int_distribution<>> rand_int;

    void outputTrajectory(std::ostream& os);
    void outputOrderParameter(std::ostream& os);
    void readxyz(const std::string filename);
    void readDiffusivity(const std::string filename);
    void runHelper(int nstep, int opt);
    void MonteCarlo();
    int CheckOverlap(int i, int j);
    void calOp();
    void calDss();
    void InitializeIndexMap();
    void InitializeEdge();
};
}

