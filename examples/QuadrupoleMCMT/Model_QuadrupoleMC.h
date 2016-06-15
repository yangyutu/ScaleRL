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
// Structures for particles
    struct coord
    {
        double x;
        double y;
    };
    struct index
    {
        int x;
        int y;
    };
    struct particle
    {
        coord center;
        double rot;
        coord edge1;
        coord edge2;
        coord edge3;
        coord edge4;
        index DisLoc;
    };
    
// Constants    
    static const int np = 300;
    static const int np3 = 900;
    static const int IndexR = 60;
    const double pi = 3.1415925025939941;
// Variables
    particle Polygon[np];
    std::vector< std::vector<int> > second;
    std::vector<std::vector<std::vector<int> > > IndexMatrix;
    int n_rows, n_cols, polygonnum, ppp;
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
    double psi6, rg, lambda, rmin, dt;
    double DiffTrans, DiffRot;
    std::uniform_real_distribution<double> rand_uniform;
    
    void outputTrajectory(std::ostream& os);
    void outputOrderParameter(std::ostream& os);
    void readxyz(const std::string filename);
    void readDiffusivity(const std::string filename);
    void runHelper(int nstep, int opt);
    void MonteCarlo();
    bool CheckOverlap(particle i, particle j);
    void calOp();
    void calDss();
    void UpdateIndex(int i);
    double Determinant(double v1, double v2, double v3, double v4);
    void UpdateEdge(int i);
};
}

