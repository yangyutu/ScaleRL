#pragma once
#include <string>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include "Model_QuadrupoleBase.h"

namespace ReinforcementLearning {  
class Model_QuadrupoleBD: public Model_QuadrupoleBase{
public:
    Model_QuadrupoleBD(){}
    Model_QuadrupoleBD(std::string filetag0, int R);
    Model_QuadrupoleBD(std::string filetag0);
    virtual ~Model_QuadrupoleBD(){ }
    virtual void run(int action);
    virtual void createInitialState();
    virtual double getRewards();
    virtual bool terminate();
protected:
    static const int np = 300;
    static const int np3 = 900;
    static const int rgdssbin = 25;
    static const int distdssbin = 50;
    const double Os_pressure = 5.8e-8;
    std::vector<std::vector<int>> nlist;
    int R, n_rows, n_cols;
    double dx1, dx2;
    int opt, nstep;
    int trajOutputInterval;
    int timeCounter,fileCounter;
    bool trajOutputFlag;
    std::ofstream trajOs, opOs;
    std::string filetag;
    double r[np3], psi6, c6, rg, lambda;
    double dssarray[rgdssbin][distdssbin];
    int dsscount[rgdssbin][distdssbin], rbin, rgbin;
    std::shared_ptr<std::uniform_int_distribution<>> rand_int;

    void outputTrajectory(std::ostream& os);
    void outputOrderParameter(std::ostream& os);
    void readxyz(const std::string filename);
    void readDiffusivity(const std::string filename);
    void runHelper(int nstep, int opt);
    void buildlist(int);

    int nxyz[np][3];
    double F[np3],D[np3];
    double randisp[np3], dsscalcu[np];
    void forces(int);
    double EMAG(double,double);
    
    void calOp();
    void calDss();
    int ecorrectflag;
    int DG;
    double a, tempr, fcm, kb, rmin,pfpp,kappa,re,rcut,fac1,fac2, dpf;
    double dssmin,dssmax, rgdsmin, delrgdsmin, distmin, deldist;
    std::shared_ptr<std::normal_distribution<>> randn_double;
};
}

