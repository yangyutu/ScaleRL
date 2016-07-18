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
#include <algorithm>

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
        std::vector<coord> edge;
        index DisLoc;
        double F,C,Rg,Psi;
    };
    
// Constants    
    static const int np = 10;   // Number of particles
    static const int IndexR = 60;   // Number of discretize states
    const double pi = 3.1415925025939941;
    
// Variables

    // Learning Map
    int n_rows, n_cols;
    double dx1, dx2;
    // Particle
    particle Polygon[np];
    double Angle, EdgeLength, a;
    double DiffTrans, DiffRot;
    int polygonnum;
    // Simulation/Control settings
    int opt;
    double nstep, dt;
    int trajOutputInterval;
    int timeCounter,fileCounter;
    std::ofstream trajOs, opOs;
    std::string filetag;
    std::ofstream file;
    std::uniform_real_distribution<double> rand_uniform;

    // OP and other state parameters
    double psi6, rg, F, lambda, C;
    double rmin, rmin2, ctestv;
    
    // Functions
    void outputTrajectory(std::ostream& os);
    void outputOrderParameter(std::ostream& os);
    void readxyz(const std::string filename);
    void runCore(int nstep, int opt);
    void MonteCarlo();
    bool CheckOverlap(particle i, particle j);
    void calPsi();
    void calRg();
    void calF();
    void calC();
};
}

