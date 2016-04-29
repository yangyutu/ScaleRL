#pragma once
#include <string>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include "Model_QuadrupoleBase.h"
#include <math.h>
#include <armadillo>

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
    const double pi = 3.1416;
    std::vector<double> RCheck, Edge;
    arma::field<std::vector<int> > IndexMap;
    int R, n_rows, n_cols, polygon, a, OverLap;
    double dx1, dx2;
    double Angle, EdgeLength;
    int nxyz[np][3],DiscretizedR[np][2];
    int opt, nstep, Index;
    int trajOutputInterval;
    int timeCounter,fileCounter;
    bool trajOutputFlag;
    std::ofstream trajOs, opOs;
    std::string filetag;
    
    double r[np3], psi6, c6, rg, lambda, rmin, dt;
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

