#include "Model_QuadrupoleMC.h"

using namespace ReinforcementLearning;

Model_QuadrupoleMC::Model_QuadrupoleMC(std::string filetag0,int Resolution, int polygon0) {
    filetag = filetag0;
    polygonnum = polygon0;
    DiffTrans = 6.362e-5;
    DiffRot = 4.772e-5;
    Angle = 2.0/polygonnum;
    rmin = 2.2;
    nstep = 10000;
    dt = 1000.0/nstep;
    stateDim = 3;
    currState.resize(stateDim);
    prevState.resize(stateDim);
    numActions = 4;
    fileCounter = 0;
    rand_normal = std::make_shared<std::normal_distribution<double>>(0.0, 1.0);
    rand_uniform = std::uniform_real_distribution<double> (0, 1);
    // 1 is the length from center to corner
    a = sin(0.5*pi-pi/polygonnum); // the length from center to edge
    EdgeLength = a*tan(pi/polygonnum); // the length of half of edge
    n_rows = Resolution;
    n_cols = Resolution;
    dx1 = 1.0/Resolution;
    dx2 = 1.0/Resolution;
}

void Model_QuadrupoleMC::createInitialState() {
    std::stringstream FileStr;
    FileStr << this->fileCounter;
//    this->readxyz("./StartMeshgridFolder/startmeshgrid" + FileStr.str() + ".txt");
    this->readxyz("./StartMeshgridFolder/startmeshgrid1.txt");
    std::stringstream ss;
    ss << this->fileCounter++;
    if (trajOs.is_open()) trajOs.close();
    if (opOs.is_open()) opOs.close();
    this->trajOs.open(filetag + "xyz_" + ss.str() + ".dat");
    this->opOs.open(filetag + "op" + ss.str() + ".dat");
    this->timeCounter = 0;
    std::vector<int> temp;
    std::vector<std::vector<int>> second;
    for (int i = 0; i < IndexR; i++){
        second.push_back(temp);
    }
    for (int i = 0; i < IndexR; i++){
        IndexMatrix.push_back(second);
    }

    for (int i = 0; i < np; i++){
        UpdateEdge(i);
        UpdateIndex(i);
    }
    for (int i = 0; i < np - 1; i++){
        for (int j = i+1; j < np; j++){
            if(this->CheckOverlap(Polygon[i],Polygon[j]) == true){
                std::cout << "Initial State Overlaps Between Particle " 
                        << i+1 << " and " << j+1<< std::endl;
                exit(2);
            }            
        }
    }
    this->runHelper(0,3);
    this->currState[0] = psi6;
    this->currState[1] = rg;
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
        linestream >> Polygon[i].center.x;
        linestream >> Polygon[i].center.y;
//        linestream >> dum;
        linestream >> Polygon[i].rot;
    }
    is.close();
}

void Model_QuadrupoleMC::UpdateEdge(int i){
        Polygon[i].edge1.x = Polygon[i].center.x + a*cos(Polygon[i].rot-pi/2);
        Polygon[i].edge1.y = Polygon[i].center.y + a*sin(Polygon[i].rot-pi/2);
        Polygon[i].edge2.x = Polygon[i].center.x + a*cos(Polygon[i].rot-pi/2+Angle*pi);
        Polygon[i].edge2.y = Polygon[i].center.y + a*sin(Polygon[i].rot-pi/2+Angle*pi);
        Polygon[i].edge3.x = Polygon[i].center.x + a*cos(Polygon[i].rot-pi/2+2*Angle*pi);
        Polygon[i].edge3.y = Polygon[i].center.y + a*sin(Polygon[i].rot-pi/2+2*Angle*pi);
}

void Model_QuadrupoleMC::UpdateIndex(int i){
    Polygon[i].DisLoc.x = floor(Polygon[i].center.x/(60/IndexR) + IndexR/2);
    Polygon[i].DisLoc.y = floor(Polygon[i].center.y/(60/IndexR) + IndexR/2);
    IndexMatrix[Polygon[i].DisLoc.x][Polygon[i].DisLoc.y].push_back(i);
}

void Model_QuadrupoleMC::outputTrajectory(std::ostream& os) {

    for (int i = 0; i < np; i++) {
        os << i << "\t";
        os << Polygon[i].center.x << "\t";
        os << Polygon[i].center.y << "\t";
        os << Polygon[i].rot<< "\t";
        os << std::endl;
    }
}
void Model_QuadrupoleMC::outputOrderParameter(std::ostream& os) {
    os << this->timeCounter << "\t";
    os << psi6 << "\t";
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
    if (currState[0] > 0.99) {return true;}
    return false;
}

void Model_QuadrupoleMC::run(int action) {
/* Takes in the control option and run the dynamics for 1s*/
    this->opt = action;
    this->outputTrajectory(this->trajOs);
    this->outputOrderParameter(this->opOs);
    this->runHelper(nstep,opt);
    this->prevState = this->currState;
    this->currState[0] = psi6;
    this->currState[1] = rg;
    this->timeCounter++;
}

void Model_QuadrupoleMC::runHelper(int nstep, int controlOpt) {
// Covert control option to lambda and run MC simulation for 1s in 1e5 steps
    switch (controlOpt)
    {
        case 0:
            lambda = 0.1;
            break;
        case 1:
            lambda = 0.2;
            break;
        case 2:
            lambda = 0.3;
            break;
        case 3:
            lambda = 1;
            break;
        default:
            lambda = 1;
            break;
    }
    for (int step = 0; step < nstep; step++) {
        this->MonteCarlo();
    }
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
    double Driftx, Drifty, RandDriftx, RandDrifty, RandRot;
// Calculate the movement of each particle (i) in 0.1ms
    for (int i = 0; i < np; i++){
// The velocity of each particle (translational and rotational)
        Driftx = -Polygon[i].center.x*lambda*dt*DiffTrans;
        Drifty = -Polygon[i].center.y*lambda*dt*DiffTrans;
        RandDriftx = (*rand_normal)(rand_generator);
        RandDrifty = (*rand_normal)(rand_generator);
        RandRot = (*rand_normal)(rand_generator);
/* The new location needs to be justified.
 * Because the location and edge information is needed to test
 * overlapping, these two data are directly replaced; other data
 * including zone information and record of zone particles are 
 * kept unchanged*/
// Store old structure in TempPolygon
        particle TempPolygon = Polygon[i];
// Update new structure in Polygon
//        Polygon[i].center.x += Driftx + RandDriftx*sqrt(DiffTrans*2.0*dt);
//        Polygon[i].center.y += Drifty + RandDrifty*sqrt(DiffTrans*2.0*dt);
        Polygon[i].center.x += Driftx;
        Polygon[i].center.y += Drifty;
        Polygon[i].rot += RandRot*sqrt(DiffRot*2.0*dt);
        UpdateEdge(i);
        Polygon[i].DisLoc.x = floor(Polygon[i].center.x/(60/IndexR) + IndexR/2);
        Polygon[i].DisLoc.y = floor(Polygon[i].center.y/(60/IndexR) + IndexR/2);
        bool OverLap = false;
        bool move = true;
//        for (int j = 0; j < np; j++){
//            if (j != i && pow((Polygon[j].DisLoc.x-Polygon[i].DisLoc.x),2.0) < 16 && pow((Polygon[j].DisLoc.y-Polygon[i].DisLoc.y),2.0) < 16 ){
//                double tempdist = sqrt(pow((Polygon[j].center.x - Polygon[i].center.x ),2.0) + 
//                       pow((Polygon[j].center.y - Polygon[i].center.y),2.0));   
//                if (tempdist < 2*a){
//                    OverLap = true;
//                }
//                if (this->CheckOverlap(Polygon[i],Polygon[j]) == true){
//                    OverLap = true;
//                }                
//            }
//        }
        for (int ii = -4; ii <= 4; ii++){
            if (Polygon[i].DisLoc.x+ii < 0 ||Polygon[i].DisLoc.x+ii >= IndexR ){
                continue;
            }
            for (int jj = -4; jj <= 4; jj++){
                if ( Polygon[i].DisLoc.y +jj < 0 ||Polygon[i].DisLoc.y >= IndexR){
                    continue;
                }
//                for (std::vector<int>::iterator kk = IndexMatrix[Polygon[i].DisLoc.x+ii][Polygon[i].DisLoc.y+jj].begin(); 
//                        kk != IndexMatrix[Polygon[i].DisLoc.x+ii][Polygon[i].DisLoc.y+jj].end(); ++kk){
                for (int kk = 0; kk < IndexMatrix[Polygon[i].DisLoc.x+ii][Polygon[i].DisLoc.y+jj].size(); kk++){
//                    int j = *kk;
                    int j = IndexMatrix[Polygon[i].DisLoc.x+ii][Polygon[i].DisLoc.y+jj].at(kk);
                    double tempdist;
                    if (j != i){
                        tempdist = sqrt(pow((Polygon[j].center.x - Polygon[i].center.x ),2.0) + 
                               pow((Polygon[j].center.y - Polygon[i].center.y),2.0));   
                        if (tempdist < 2*a){
                            OverLap = true;
                        }
                        if (this->CheckOverlap(Polygon[i],Polygon[j]) == true){
                            OverLap = true;
                        }
                    }
                }
            }
        }
// If no overlap, accept the move at the ratio of Boltzmann distribution        
       double PotentialDiff = (pow(Polygon[i].center.x,2.0) + pow(Polygon[i].center.y,2.0)) 
       - (pow(TempPolygon.center.x,2.0) + pow(TempPolygon.center.y,2.0));
       if (PotentialDiff > 0){
           double RanGen = rand_uniform(rand_generator);
           double BoltzmannDist = exp(-0.5*lambda*
           ((pow(Polygon[i].center.x,2.0) + pow(Polygon[i].center.y,2.0)) - 
           (pow(TempPolygon.center.x,2.0) + pow(TempPolygon.center.y,2.0))));
           if (RanGen >= BoltzmannDist){
               move = false;
           }
       }
// If no overlap, update IndexMap
        if (OverLap == false && move == true){
            IndexMatrix[Polygon[i].DisLoc.x][Polygon[i].DisLoc.y].push_back(i);
            IndexMatrix[TempPolygon.DisLoc.x][TempPolygon.DisLoc.y].erase(IndexMatrix[TempPolygon.DisLoc.x][TempPolygon.DisLoc.y].begin());            
        } else {
            Polygon[i] = TempPolygon;
        }
    }
}

bool Model_QuadrupoleMC::CheckOverlap(particle i, particle j){
/* check overlap by filling in the indices of two particles needed to be checked;
 * the result is 0 if not overlap and 10 otherwise. The checking requires r matrix and 
 * edge matrix information */
    double DiffEdge1, DiffEdge2, Dist, Det, Det1, Det2, Frac1, Frac2;
    double Edge1x[3] = {i.edge1.x,i.edge2.x,i.edge3.x};
    double Edge1y[3] = {i.edge1.y,i.edge2.y,i.edge3.y};
    double Edge2x[3] = {j.edge1.x,j.edge2.x,j.edge3.x};
    double Edge2y[3] = {j.edge1.y,j.edge2.y,j.edge3.y};
    for (int ii = 0; ii < polygonnum; ii++){
        for (int jj = 0; jj < polygonnum; jj++){
            DiffEdge1 = Edge2x[jj] - Edge1x[ii];
            DiffEdge2 = Edge2y[jj] - Edge1y[ii];
            Dist = sqrt(DiffEdge1*DiffEdge1 + DiffEdge2*DiffEdge2);
            if (Dist < 2*EdgeLength){
                Det = -cos(i.rot+ii*Angle*pi)*sin(j.rot+jj*Angle*pi)
                        + cos(j.rot+jj*Angle*pi)*sin(i.rot+ii*Angle*pi);
                if (Det != 0.0){
                    Det1 = -DiffEdge1*sin(j.rot+jj*Angle*pi) + DiffEdge2*cos(j.rot+jj*Angle*pi);
                    Det2 = DiffEdge2*cos(i.rot+ii*Angle*pi) - DiffEdge1*sin(i.rot+ii*Angle*pi); 
                    Frac1 = fabs(Det1/Det);
                    Frac2 = fabs(Det2/Det);
//                    file << i+1 << "\t" << j+1 <<"\t" << ii+1 << "\t" << jj+1 << "\t"
//                            << DiffEdge1 << "\t" << DiffEdge2 << "\t" << Dist << "\t"
//                            << Det << "\t" << Det1 << "\t" << Det2 << "\t"
//                            << Frac1 << "\t" << Frac2 << std::endl;
                    if (Frac1 <= EdgeLength && Frac2 <= EdgeLength){
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

void Model_QuadrupoleMC::calOp() {

    int nb[np];
    double rx[np], ry[np], phi[np];
    double rxij, ryij, psir[np], psii[np], scale;
    double rgmean, xmean, ymean, accumpsi6r, accumpsi6i;
    
    if (polygonnum == 3) {scale = 6.0;} 
    else if (polygonnum == 4){scale = 4.0;} 
    else if (polygonnum == 6){scale = 6.0;} 
    else {scale = 1.0;}
    
    for (int i = 0; i < np; i++) {
        rx[i] = Polygon[i].center.x;
        ry[i] = Polygon[i].center.y;
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
                if (RP <= rmin) {
                    nb[i] += 1;
                    psir[i] += cos(scale * Polygon[j].rot);
                    psii[i] += sin(scale * Polygon[j].rot);
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
        rgmean = rgmean + pow((rx[i] - xmean),2) + pow((rx[i] - xmean),2);
        rgmean = rgmean + pow((ry[i] - ymean),2) + pow((ry[i] - xmean),2);
    }
    rgmean /= np;
    rg = sqrt(rgmean);
    rg = 1.0 - (rg - 6.0)/24.0;
}