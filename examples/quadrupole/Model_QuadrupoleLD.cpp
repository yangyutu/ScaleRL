#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h>

#include "Model_QuadrupoleLD.h"
using namespace ReinforcementLearning;

Model_QuadrupoleLD::Model_QuadrupoleLD(std::string filetag0) {
    filetag = filetag0;
    currState.resize(2);
    prevState.resize(2);
    xdelop = 0.025;
    ydelop = -50.0;
    xbinNum = 40;
    ybinNum = 156;
    xminop = 0.0;
    yminop = 24000.0;
    dt = 0.10; //s
    stateDim = 2;
    numActions = 3;
    trajOutputInterval = 10;
    fileCounter = 0;
    std::array<Array4D_type::index, 4> dims1 = {numActions, 2, xbinNum, ybinNum};
    std::array<Array2D_type::index, 4> dims2 = {xbinNum, ybinNum};
    v = std::make_shared<Array4D_type>(dims1);
    d = std::make_shared<Array4D_type>(dims1);
    hasValue = std::make_shared<Array2D_type>(dims2);
    rand_normal = std::make_shared<std::normal_distribution<double>>(0.0, 1.0);


    std::string filenames[3];
    // the 0.3v
    filenames[0] = "templatevandd1.txt";
    // the 0.45v	
    filenames[1] = "templatevandd3.txt";
    // the 0.7v
    filenames[2] = "templatevandd5.txt";
    load_V_D(filenames);
}

void Model_QuadrupoleLD::load_V_D(std::string filenames[]) {


    
    double vxtemp, vytemp, dxxtemp, dyytemp;
    int hasValuetemp;
    double dum;

    for (int i = 0; i < xbinNum; i++) {
        for (int j = 0; j < ybinNum; j++) {
            (*hasValue)[i][j] = 1;
        }
    }
    for (int m = 0; m < this->numActions; m++) {
        std::ifstream is;
        std::string line;
        is.open(filenames[m]);
        for (int i = 0; i < xbinNum; i++) {
            for (int j = 0; j < ybinNum; j++) {
                getline(is, line);
                std::stringstream linestream(line);
                linestream >> dum;
                linestream >> dum;
                linestream >> dum;
                linestream >> vxtemp;
                linestream >> vytemp;
                linestream >> dxxtemp;
                linestream >> dum;
                linestream >> dum;
                linestream >> dyytemp;
                linestream >> hasValuetemp;
                (*v)[m][0][i][j] = vxtemp * 1000.0;
                (*v)[m][1][i][j] = vytemp * 1000.0;
                (*d)[m][0][i][j] = dxxtemp * 1000.0;
                (*d)[m][1][i][j] = dyytemp * 1000.0;
                if (!hasValuetemp) (*hasValue)[i][j] = 0;
            }
        }
    }

    std::ofstream os;
    os.open("stateSpace.txt");
    int count = 0;
    for (int i = 0; i < xbinNum; i++) {
        for (int j = 0; j < ybinNum; j++) {
            if ((*hasValue)[i][j]) {
                count++;
                os << i << "\t" << j << "\t";
                os << (xminop + xdelop*(i+0.5)) << "\t";
                os << (yminop + ydelop*(j+0.5)) << std::endl;
                
            }
        }
    }
    os.close();

    std::cout << "there are total states of " << count << std::endl;
}

void Model_QuadrupoleLD::run(int action) {

    if (this->timeCounter == 0 || (this->timeCounter+1)%trajOutputInterval == 0){
        this->outputTrajectory(this->trajOs, action);
    }
    int m = (int) ((currState[0] - xminop) / xdelop) + 1;
    int n = (int) ((currState[1] - yminop) / ydelop) + 1;
    if (m < 0) m = 0;
    if (n < 0) n = 0;
    if (m >= xbinNum) m = xbinNum - 1;
    if (n >= ybinNum) n = ybinNum - 1;
    double vxtemp = (*v)[action][0][m][n];
    double vytemp = (*v)[action][1][m][n];
    double dxxtemp = (*d)[action][0][m][n];
    double dyytemp = (*d)[action][1][m][n];
    double tempx = currState[0] + vxtemp * dt +  ((*rand_normal)(rand_generator)) * sqrt(2 * dxxtemp *dt);
    double tempy = currState[1] + vytemp * dt +  ((*rand_normal)(rand_generator)) * sqrt(2 * dyytemp *dt);
    
    if (tempx < 0.0) tempx = -tempx;
    
    m = (int) ((currState[0] - xminop) / xdelop) + 1;
    n = (int) ((currState[1] - yminop) / ydelop) + 1;
    if (m < 0) m = 0;
    if (n < 0) n = 0;
    if (m >= xbinNum) m = xbinNum - 1;
    if (n >= ybinNum) n = ybinNum - 1;
    
    if ((*hasValue)[m][n]) {
        currState[0] = tempx;
        currState[1] = tempy; 

    } else {
        tempx *=0.9;
        tempy *=1.0;
        currState[0] = tempx;
        currState[1] = tempy;        
        
        std::cerr << "in forbidden region" << std::endl;    
    }

   this->timeCounter++;    

}

void Model_QuadrupoleLD::createInitialState() {
    this->currState[0] = 0.05;
    this->currState[1] = 20000;
    std::stringstream ss;
    ss << this->fileCounter++;
    
    this->trajOs.open(filetag+ss.str()+".dat");
    this->timeCounter = 0;
    
}

void Model_QuadrupoleLD::outputTrajectory(std::ostream& os, int action){
   os << this->timeCounter * dt << "\t";
    for (int s = 0; s < this->stateDim; s++) {
        os << this->currState[s] << "\t";
    }
    os << action << std::endl;
}

double Model_QuadrupoleLD::getRewards() {
    if (this->terminate()) {
        this->reward = 0;
        return reward;
    } else {
        this->reward = -(1- currState[0]);
        return reward;
    }
}

bool Model_QuadrupoleLD::terminate() {
    // the model will stop if \psi6 > 0.95, which is perfect crystal
    if (currState[0] > 0.90){
        this->trajOs.close();
        return true;
    }
    return false;
}
