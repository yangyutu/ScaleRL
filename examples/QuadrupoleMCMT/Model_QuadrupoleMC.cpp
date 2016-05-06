#include "Model_QuadrupoleMC.h"

using namespace ReinforcementLearning;

Model_QuadrupoleMC::Model_QuadrupoleMC(std::string filetag0,int R, int polygon0) {
    filetag = filetag0;
    polygon = polygon0;
    DiffTrans = 0.00006362;
    DiffRot = 0.00004772;
    Angle = 2.0/polygon;
    EdgeLength = tan(pi/polygon);
    rmin = 2.2;
    nstep = 10000;
    dt = 1000.0/nstep;
    stateDim = 3;
    currState.resize(stateDim);
    prevState.resize(stateDim);
    numActions = 4;
    fileCounter = 0;
    rand_normal = std::make_shared<std::normal_distribution<double>>(0.0, 1.0);
    a = sin(0.5*pi-pi/polygon);
    n_rows = R;
    n_cols = R;
    dx1 = 1.0/R;
    dx2 = 1.0/R;
}


void Model_QuadrupoleMC::run(int action) {
/* Takes in the control option and run the dynamics for 1s*/
    this->opt = action;
    this->outputTrajectory(this->trajOs);
    this->outputOrderParameter(this->opOs);
    this->runHelper(nstep,opt);
    c6 = c6 / 5.6;
    this->prevState = this->currState;
    this->currState[0] = psi6;
    this->currState[1] = c6;
    this->timeCounter++;
}

void Model_QuadrupoleMC::createInitialState() {
    std::stringstream FileStr;
    FileStr << this->fileCounter;
//    this->readxyz("./StartMeshgridFolder/startmeshgrid" + FileStr.str() + ".txt");
    this->readxyz("./StartMeshgridFolder/startmeshgrid1.txt");
    std::stringstream ss;
    std::cout << "model initialize at round " << fileCounter << std::endl;
    ss << this->fileCounter++;
    if (trajOs.is_open()) trajOs.close();
    if (opOs.is_open()) opOs.close();
    this->trajOs.open(filetag + "xyz_" + ss.str() + ".dat");
    this->opOs.open(filetag + "op" + ss.str() + ".dat");
    this->timeCounter = 0;
    
    
    this->InitializeEdge();

    for (int i = 0; i < np - 1; i++){
        for (int j = i+1; j < np; j++){
            if(this->CheckOverlap(i,j) > 2){
                std::cout << "Initial State Overlaps Between Particle " 
                        << i+1 << " and " << j+1<< std::endl;
                exit(2);
            }            
        }
    }
    std::cout << " No Overlapping at Starting State" << std::endl;
    this->runHelper(0,3);
    this->InitializeIndexMap();
    c6 = c6 / 5.6;
    this->currState[0] = psi6;
    this->currState[1] = c6;
}
void Model_QuadrupoleMC::InitializeIndexMap(){
    IndexMap.set_size(IndexR,IndexR);
    for (int i = 0; i < np; i++){
        DiscretizedR[i][0] = ceil(r[i*3+0]/(60/IndexR)) + IndexR/2;
        DiscretizedR[i][1] = ceil(r[i*3+1]/(60/IndexR)) + IndexR/2;
        IndexMap(DiscretizedR[i][0],DiscretizedR[i][1]).push_back(i);
    }
}

void Model_QuadrupoleMC::InitializeEdge(){
    for (int i = 0; i < np; i++){
        for (int j = 0; j < polygon; j++){
            Edge.push_back(r[i*3] + a*cos(r[i*3+2] - 0.5*pi + j*Angle*pi));
            Edge.push_back(r[i*3+1] + a*sin(r[i*3+2] - 0.5*pi + j*Angle*pi));
        }
    }
}

void Model_QuadrupoleMC::outputTrajectory(std::ostream& os) {

    for (int i = 0; i < np; i++) {
        os << i << "\t";
        os << r[3 * i] << "\t";
        os << r[3 * i + 1] << "\t";
        os << r[3 * i + 2]<< "\t";
        os << std::endl;
    }
}

void Model_QuadrupoleMC::outputOrderParameter(std::ostream& os) {
    os << this->timeCounter << "\t";
    os << psi6 << "\t";
    os << c6 << "\t";
    os << rg << "\t";
    os << opt << "\t";
    os << lambda << "\t";
    os << std::endl;
}

double Model_QuadrupoleMC::getRewards() {
    if (this->terminate()) {
        this->reward = 0;
        return reward;
    } else {
        this->reward = -(1 - currState[0]);
        return reward;
    }
}

bool Model_QuadrupoleMC::terminate() {
    if (currState[0] > 0.90) {return true;}
    return false;
}

void Model_QuadrupoleMC::readxyz(const std::string filename) {
    std::ifstream is;
    is.open(filename.c_str());
    std::string line;
    double dum;
    for (int i = 0; i < np; i++) {
        getline(is, line);
        std::stringstream linestream(line);
        linestream >> dum;
        linestream >> r[3 * i];
        linestream >> r[3 * i + 1];
        linestream >> dum;
        linestream >> r[3 * i + 2];
    }
    is.close();
}

void Model_QuadrupoleMC::runHelper(int nstep, int controlOpt) {
// Covert control option to lambda and run MC simulation for 1s in 1e5 steps
    if (controlOpt == 0) {lambda = 0.1;}
    else if (controlOpt == 1) {lambda = 0.2;}
    else if (controlOpt == 2) {lambda = 0.3;}
    else {lambda = 1;}
//    file.open("Traj.dat");
    for (int step = 0; step < nstep; step++) {
        this->MonteCarlo();
    }
//    file.close();
    this->calOp();
}

void Model_QuadrupoleMC::MonteCarlo(){
/* Use MC method to simulate the trajectory of particles for 1s.
 * The idea is to update the location and motion of each particle 
 * ever 0.01ms.
 * For each 0.01ms, each particle is calculated its movement,
 * if the motion is allowed (new location not overlapping with
 * other particles, then the movement is kept; otherwise, this
 * movement is discarded, and this particle is considered to stay
 * still for this 0.1ms*/
    int DiscretizedRNew[2],Index;
    double Driftx, Drifty, RandDriftx, RandDrifty, RandRot, TempR[3];
    std::vector<double> TempEdge;
// Calculate the movement of each particle (i) in 0.1ms
    for (int i = 0; i < np; i++){
//        file << i+1 << "\t" << r[i*3+0] << "\t" << r[i*3+1] << std::endl;
// The velocity of each particle (translational and rotational)
        Driftx = -r[3*i]*lambda*dt*DiffTrans;
        Drifty = -r[3*i+1]*lambda*dt*DiffTrans;
        RandDriftx = (*rand_normal)(rand_generator);
        RandDrifty = (*rand_normal)(rand_generator);
        RandRot = (*rand_normal)(rand_generator);
/* The new location needs to be justified.
 * Because the location and edge information is needed to test
 * overlapping, these two data are directly replaced; other data
 * including zone information and record of zone particles are 
 * kept unchanged*/
// Record the old location of particle i as TempR
        TempR[0] = r[i*3+0];
        TempR[1] = r[i*3+1];
        TempR[2] = r[3*i+2];
// Update the new location into location matrix r
        r[i*3+0] = r[i*3+0] + Driftx + RandDriftx*sqrt(DiffTrans*2.0*dt);
        r[i*3+1] = r[i*3+1] + Drifty + RandDrifty*sqrt(DiffTrans*2.0*dt);
        r[i*3+2] = r[i*3+2] + RandRot*sqrt(DiffRot*2.0*dt);
//        r[i*3+0] = r[i*3+0] + Driftx;
//        r[i*3+1] = r[i*3+1] + Drifty;
//        r[i*3+2] = r[i*3+2];
        
// Record the old edge as TempEdge
        for (int kk = 0; kk < polygon; kk++){
            TempEdge.push_back(Edge.at(i*2*polygon+kk*2+0));
            TempEdge.push_back(Edge.at(i*2*polygon+kk*2+1));                
        }
// Update the new edge into edge matrix
        for (int kk = 0; kk < polygon; kk++){
            Edge.at(i*2*polygon+kk*2+0) = r[i*3+0]+a*cos(r[i*3+2]-0.5*pi+(kk)*Angle*pi);
            Edge.at(i*2*polygon+kk*2+1) = r[i*3+1]+a*sin(r[i*3+2]-0.5*pi+(kk)*Angle*pi);                
        }        
// Record new zone info into temp DiscretizedRNew
        DiscretizedRNew[0] = ceil(r[i*3+0]/(60/IndexR)) + IndexR/2;
        DiscretizedRNew[1] = ceil(r[i*3+1]/(60/IndexR)) + IndexR/2;
/* Test overlapping of particles that are in the 8 zones around
 * the location of new particle i location */
        int OverLapTot = 0;
//        for (int ii = -1; ii <= 1; ii++){
//            for (int jj = -1; jj <= 1; jj++){
//                for (int kk = 0; kk < IndexMap(DiscretizedRNew[0]+ii,DiscretizedRNew[1]+jj).size();kk++){
//                    Index = IndexMap(DiscretizedRNew[0]+ii,DiscretizedRNew[1]+jj).at(kk);
//                    if (Index != i && (this->CheckOverlap(i,Index) > 1)){
//                        OverLapTot = 10;
//                    }                    
//                }
//            }
//        }
        for (int ii = 0; ii < np; ii++){
            if (ii != i && (this->CheckOverlap(i,ii) > 1)){
                OverLapTot = 10;
                break;
            }
        }
// If no overlap, this movement is allowed; update IndexMap and DiscretizedR 
        if (OverLapTot < 1){
            IndexMap(DiscretizedRNew[0],DiscretizedRNew[1]).push_back(IndexMap(DiscretizedR[i][0],DiscretizedR[i][1]).at(0));
            IndexMap(DiscretizedR[i][0],DiscretizedR[i][1]).erase(IndexMap(DiscretizedR[i][0],DiscretizedR[i][1]).begin());
            DiscretizedR[i][0] = DiscretizedRNew[0];
            DiscretizedR[i][1] = DiscretizedRNew[1];
        }
// If overlap, this movement is not allowed; restore r matrix and Edge matrix
        else{
            r[i*3+0] = TempR[0];
            r[i*3+1] = TempR[1];
            r[i*3+2] = TempR[2];
            for (int kk = 0; kk < polygon; kk++){
                Edge.at(i*2*polygon+kk*2+0) = TempEdge.at(kk*2+0);
                Edge.at(i*2*polygon+kk*2+1) = TempEdge.at(kk*2+1);
            }
        }
        TempEdge.clear();
    }
}

int Model_QuadrupoleMC::CheckOverlap(int i, int j){
/* check overlap by filling in the indices of two particles needed to be checked;
 * the result is 0 if not overlap and 10 otherwise. The checking requires r matrix and 
 * edge matrix information */
    double DiffEdge1, DiffEdge2, Dist, Det, Det1, Det2, Frac1, Frac2;
    for (int ii = 0; ii < polygon; ii++){
        for (int jj = 0; jj < polygon; jj++){
            DiffEdge1 = Edge.at(j*2*polygon+jj*2+0) - Edge.at(i*2*polygon+ii*2+0);
            DiffEdge2 = Edge.at(j*2*polygon+jj*2+1) - Edge.at(i*2*polygon+ii*2+1);
            Dist = sqrt(DiffEdge1*DiffEdge1 + DiffEdge2*DiffEdge2);
/* First check if the distance between particles are too close; if so no other criteria
 * are needed; otherwise edges needed to be checked */
            if (Dist < 2*EdgeLength){
                Det = -cos(r[i*3+2]+ii*Angle*pi)*sin(r[j*3+2]+jj*Angle*pi)
                        + cos(r[j*3+2]+jj*Angle*pi)*sin(r[i*3+2]+ii*Angle*pi);
                if (Det != 0.0){
                    Det1 = -DiffEdge1*sin(r[j*3+2]+jj*Angle*pi) + DiffEdge2*cos(r[j*3+2]+jj*Angle*pi);
                    Det2 = DiffEdge2*cos(r[i*3+2]+ii*Angle*pi) - DiffEdge1*sin(r[i*3+2]+ii*Angle*pi); 
                    Frac1 = fabs(Det1/Det);
                    Frac2 = fabs(Det2/Det);
//                    file << i+1 << "\t" << j+1 <<"\t" << ii+1 << "\t" << jj+1 << "\t"
//                            << DiffEdge1 << "\t" << DiffEdge2 << "\t" << Dist << "\t"
//                            << Det << "\t" << Det1 << "\t" << Det2 << "\t"
//                            << Frac1 << "\t" << Frac2 << std::endl;
                    if (Frac1 <= EdgeLength && Frac2 <= EdgeLength){
                        return 10;
                    }
                }
            }
        }
    }
    return 0;
}

void Model_QuadrupoleMC::calOp() {

    int nb[np], con[np];
    double rx[np], ry[np];
    double rxij, ryij, theta, psir[np], psii[np], numer, denom, testv, ctestv;
    double rgmean, xmean, ymean, accumpsi6r, accumpsi6i;
    ctestv = 0.32;
    for (int i = 0; i < np; i++) {
        rx[i] = r[i*3+0];
        ry[i] = r[i*3+1];
    }
    for (int i = 0; i < np; i++) {
        nb[i] = 0;
        psir[i] = 0.0;
        psii[i] = 0.0;
        for (int j = 0; j < np; j++) {
            if (i != j) {
                rxij = rx[j] - rx[i];
                ryij = ry[j] - ry[i];
                double RP = sqrt(rxij * rxij + ryij * ryij);
                if (RP < rmin) {
                    nb[i] += 1;
                    theta = std::atan2(ryij, rxij);
                    psir[i] += cos(6 * theta);
                    psii[i] += sin(6 * theta);
                }
            }        
        } 
            if (nb[i] > 0) {
                psir[i] /=  nb[i];
                psii[i] /=  nb[i];
            }
    }
    psi6 = 0;
    accumpsi6r = 0;
    accumpsi6i = 0;
    for (int i = 0; i < np; i++) {
        accumpsi6r = accumpsi6r + psir[i];
        accumpsi6i = accumpsi6i + psii[i];
    }
    accumpsi6r = accumpsi6r / np;
    accumpsi6i = accumpsi6i / np;
    psi6 = sqrt(accumpsi6r * accumpsi6r + accumpsi6i * accumpsi6i);
    c6 = 0.0;
    for (int i = 0; i < np; i++) {
        con[i] = 0;
        for (int j = 0; j < np; j++) {
            rxij = rx[j] - rx[i];
            ryij = ry[j] - ry[i];
            double rp = sqrt(rxij * rxij + ryij * ryij);
            if ((i != j)&&(rp <= rmin)) {

                numer = psir[i] * psir[j] + psii[i] * psii[j];
                double temp = psii[i] * psir[j] - psii[j] * psir[i];
                denom = sqrt(numer * numer + temp*temp);
                testv = numer / denom;
                if (testv >= ctestv) {
                    con[i] += 1;
                }
            }
        }
        c6 = c6 + con[i];
    }
    c6 /= np;
    //      calculate Rg
    xmean = 0;
    ymean = 0;
    for (int i = 0; i < np; i++) {
        xmean = xmean + rx[i];
        ymean = ymean + ry[i];
    }
    xmean /= np;
    ymean /= np;
    rgmean = 0;
    for (int i = 0; i < np; i++) {
        rgmean = rgmean + (rx[i] - xmean)*(rx[i] - xmean);
        rgmean = rgmean + (ry[i] - ymean)*(ry[i] - xmean);
    }
    rgmean /= np;
    rgmean = sqrt(rgmean);
    rg = rgmean;
}
