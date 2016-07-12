#include "Model_QuadrupoleMC.h"

using namespace ReinforcementLearning;

Model_QuadrupoleMC::Model_QuadrupoleMC(std::string filetag0,int Resolution, int polygon0) {
    filetag = filetag0;
    polygonnum = polygon0;
    DiffTrans = 6.362e-5;
    DiffRot = 4.772e-5;
    Angle = 2.0/polygonnum;
    rmin = 2.2;
    rmin2 = 1.5;
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
    a = sin(0.5*pi-pi/polygonnum); // the length from center to middle of edge
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
    this->readxyz("./StartMeshgridFolder/startmeshgrid1.txt"); // load initial particle configuration
    std::stringstream ss;
    ss << this->fileCounter++;
    if (trajOs.is_open()) trajOs.close();
    if (opOs.is_open()) opOs.close();
    this->trajOs.open(filetag + "xyz_" + ss.str() + ".dat"); // file for coordinate record
    this->opOs.open(filetag + "op" + ss.str() + ".dat"); // file for order parameter record
    this->timeCounter = 0;
    coord temp;
    for (int i = 0; i < np; i++){ // Coordinate of edges and particle zone
        Polygon[i].DisLoc.x = floor(Polygon[i].center.x/(60/IndexR) + IndexR/2);
        Polygon[i].DisLoc.y = floor(Polygon[i].center.y/(60/IndexR) + IndexR/2);
        for (int j = 0; j < 3; j++){
            temp.x = Polygon[i].center.x + a*cos(Polygon[i].rot - pi/2 + j*Angle*pi);
            temp.y = Polygon[i].center.y + a*sin(Polygon[i].rot - pi/2 + j*Angle*pi);
            Polygon[i].edge.push_back(temp);
        }
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
    this->currState[2] = F;
}

void Model_QuadrupoleMC::readxyz(const std::string filename) {
    // Load initial configuration as coordinate of polygon center and the rotation of one edge
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
        linestream >> dum;
        linestream >> Polygon[i].rot;
    }
    is.close();
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
    os << F << "\t";
    os << opt << "\t";
    os << lambda << "\t";
    os << std::endl;
}

double Model_QuadrupoleMC::getRewards() {
    if (this->terminate()) {
        this->reward = 0;
        return reward;
    } else {
        this->reward = -(1 - (currState[0]+currState[2])/2);
        return reward;
    }
}

bool Model_QuadrupoleMC::terminate() {
    if (currState[0] > 0.99 && currState[2] > 0.99) {return true;}
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
    this->currState[2] = F;
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
/* Simulation the motions of all particles in one second by 1e6 sweeps. In each sweep, 
 * one particle is moved at a time.
 * 
 * The movement of a particle i is immediately updated in Polygon[i], while the original 
 * i is stored in TempPolygon. 
 * 
 * Polygon[i] will be reset to TempPolygon when:
 *  1.  the movement causes overlap between particle i and other particles, or
 *  2.  the movement causes an increase in the total energy of the ensemble, and the
 *      movement is rejected following a Boltzmann distribution. */
    
    double Driftx, Drifty, RandDriftx, RandDrifty, RandRot;
    for (int i = 0; i < np; i++){
        // Save the current position into TempPolygon, and calculate the new position
        particle TempPolygon = Polygon[i];
        Driftx = -Polygon[i].center.x*lambda*dt*DiffTrans;
        Drifty = -Polygon[i].center.y*lambda*dt*DiffTrans;
        RandDriftx = (*rand_normal)(rand_generator);
        RandDrifty = (*rand_normal)(rand_generator);
        RandRot = (*rand_normal)(rand_generator);
        Polygon[i].center.x = Polygon[i].center.x + Driftx ;
        Polygon[i].center.y = Polygon[i].center.y + Drifty ;
//        Polygon[i].rot = Polygon[i].rot + RandRot*sqrt(DiffRot*2.0*dt);
        Polygon[i].DisLoc.x = floor(Polygon[i].center.x/(60/IndexR) + IndexR/2);
        Polygon[i].DisLoc.y = floor(Polygon[i].center.y/(60/IndexR) + IndexR/2);
        for (int edge = 0; edge < 3; edge++){
            Polygon[i].edge.at(edge).x = Polygon[i].center.x + a*cos(Polygon[i].rot - pi/2 + edge*Angle*pi);
            Polygon[i].edge.at(edge).y = Polygon[i].center.y + a*sin(Polygon[i].rot - pi/2 + edge*Angle*pi);
        }
        // Check overlap, return false if no overlap
        bool OverLap = false;
        for (int ii = 0; ii < np; ii++){
            if (ii != i && fabs(Polygon[ii].DisLoc.x -Polygon[i].DisLoc.x) <= 3  
                    && fabs(Polygon[ii].DisLoc.y -Polygon[i].DisLoc.y) <= 3){
                double tempdist = sqrt(pow((Polygon[ii].center.x - Polygon[i].center.x ),2.0) + 
                       pow((Polygon[ii].center.y - Polygon[i].center.y),2.0));   
                if (tempdist < 2*a){
                    OverLap = true;
                    break;
                }else if (this->CheckOverlap(Polygon[i],Polygon[ii]) == true){
                    OverLap = true;
                    break;
                }                
            }
            if (OverLap == true){break;}
        }
        // If no overlap, test potential energy change
        if (OverLap == false){
            double PotentialDiff = 
            (pow(Polygon[i].center.x,2.0) + pow(Polygon[i].center.y,2.0)) 
           - (pow(TempPolygon.center.x,2.0) + pow(TempPolygon.center.y,2.0));
            // If the potential decreases the movement is accepted; otherwise
            // the step is accepted at a probability following Boltzmann distribution
            if (PotentialDiff > 0){
               double RanGen = rand_uniform(rand_generator);
               double BoltzmannDist = 
               exp(-0.5*lambda*((pow(Polygon[i].center.x,2.0) + pow(Polygon[i].center.y,2.0))
               - (pow(TempPolygon.center.x,2.0) + pow(TempPolygon.center.y,2.0))));
               if (RanGen < BoltzmannDist){
                   Polygon[i] = TempPolygon;
               }
            }
        } else {
            Polygon[i] = TempPolygon;
        }
    }
}

bool Model_QuadrupoleMC::CheckOverlap(particle i, particle j){
// check overlap between particle i and j, return false if no overlap
    double Diffx, Diffy, Dist, Det, Det1, Det2, Frac1, Frac2;
    for (int ii = 0; ii < polygonnum; ii++){
        for (int jj = 0; jj < polygonnum; jj++){
            Diffx = j.edge.at(jj).x - i.edge.at(ii).x;
            Diffy = j.edge.at(jj).y - i.edge.at(ii).y;
            Dist = sqrt(Diffx*Diffx + Diffy*Diffy);
            if (Dist < 2*EdgeLength){
                Det = -cos(i.rot+ii*Angle*pi)*sin(j.rot+jj*Angle*pi)
                        + cos(j.rot+jj*Angle*pi)*sin(i.rot+ii*Angle*pi);
                if (Det != 0.0){
                    Det1 = -Diffx*sin(j.rot + jj*Angle*pi) + Diffy*cos(j.rot + jj*Angle*pi);
                    Det2 = Diffy*cos(i.rot + ii*Angle*pi) - Diffx*sin(i.rot + ii*Angle*pi); 
                    Frac1 = fabs(Det1/Det);
                    Frac2 = fabs(Det2/Det);
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

    int neighbor[np];
    double psir[np], psii[np], scale;
    double rgmean, accumpsi6r, accumpsi6i;
    
    if (polygonnum == 3) {scale = 6.0;} 
    else if (polygonnum == 4){scale = 4.0;} 
    else if (polygonnum == 6){scale = 6.0;} 
    else {scale = 1.0;}
    // Calculate Psi6
    for (int i = 0; i < np; i++){
        neighbor[i] = 0;
        psir[i] = 0.0;
        psii[i] = 0.0;
        accumpsi6r = 0;
        accumpsi6i = 0;
    }
    for (int i = 0; i < np-1; i++) {
        for (int j = i+1; j < np; j++) {
            double RP = sqrt(pow((Polygon[i].center.x-Polygon[j].center.x),2.0) 
            + pow((Polygon[i].center.y-Polygon[j].center.y),2.0));
            if (RP <= rmin) {
                neighbor[i] += 1;
                neighbor[j] += 1;
                psir[i] += cos(scale * Polygon[j].rot);
                psii[i] += sin(scale * Polygon[j].rot);
                psir[j] += cos(scale * Polygon[i].rot);
                psii[j] += sin(scale * Polygon[i].rot);
            }
        }
    }
    for (int i = 0; i < np; i++){
        if (neighbor[i] != 0) {
            psir[i] /=  neighbor[i];
            psii[i] /=  neighbor[i];
        }
    }

    for (int i = 0; i < np; i++) {
        accumpsi6r += psir[i];
        accumpsi6i += psii[i];
    }
    accumpsi6r /= np;
    accumpsi6i /= np;
    psi6 = sqrt(pow(accumpsi6r,2.0) + pow(accumpsi6i,2.0));
    
    // Calculate Rg
    coord mean;
    mean.x = 0;
    mean.y = 0;
    rgmean = 0;
    for (int i = 0; i < np; i++) {
        mean.x += Polygon[i].center.x;
        mean.y += Polygon[i].center.x;
    }
    mean.x /= np;
    mean.y /= np;
    for (int i = 0; i < np; i++) {
        rgmean += pow((Polygon[i].center.x - mean.x),2) + pow((Polygon[i].center.y - mean.y),2);
    }
    rgmean /= np;
    rg = sqrt(rgmean);
    
    /* New Order Parameter F.
     * 1.   The center-center distance must be shorter than 1.5
     * 2.   The angle between center-center and edge-norm vector is smaller than 0.7559
     * 
     */
    coord VecCC, Norm1, VecEE, Norm2, VecE1, VecE2;
    double DistCC,angle, Elevation, Misalign, Overlap;
//    Initialize LocF: F values of each particles, and NeighborF: number of surroudings
    F = 0;
    for (int i = 0; i < np; i++){
        LocF[i] = 0.0;
        NeighborF[i] = 0;        
    }
    for (int i = 0; i < np; i++){
        for (int j = i+1; j < np; j++){
            // Center-center distance of two particles, if close enough, calculate overlap
            DistCC = sqrt(pow((Polygon[i].center.x-Polygon[j].center.x),2.0) 
                    + pow((Polygon[i].center.y-Polygon[j].center.y),2.0));
            if (DistCC < 1.5){
//                Vector connecting two particle centers
                VecCC.x = (Polygon[j].center.x-Polygon[i].center.x);
                VecCC.y = (Polygon[j].center.y-Polygon[i].center.y);
                for (int ii = 0; ii < 3; ii++){ // each edge of particle 1
//                    normal vector of each edge of particle 1
                    Norm1.x = (Polygon[i].edge.at(ii).x - Polygon[i].center.x)/a;
                    Norm1.y = (Polygon[i].edge.at(ii).y - Polygon[i].center.y)/a;
//                    cosine of angle between Norm1 and VecCC
                    angle = (Norm1.x*VecCC.x + Norm1.y*VecCC.y)
                            /sqrt(VecCC.x*VecCC.x + VecCC.y*VecCC.y);
                    if (angle >= 0.7559){
//                        theta = Polygon[i].rot + ii*Angle*pi;
                        for (int jj = 0; jj < 3; jj++){ // each edge of particle 2
                            VecEE.x = Polygon[j].edge.at(jj).x - Polygon[i].edge.at(ii).x;
                            VecEE.y = Polygon[j].edge.at(jj).y - Polygon[i].edge.at(ii).y;
                            Norm2.x = (Polygon[j].edge.at(jj).x - Polygon[j].center.x)/a;
                            Norm2.y = (Polygon[j].edge.at(jj).y - Polygon[j].center.y)/a;
                            Elevation = std::max(
                                    fabs(VecEE.x*Norm1.x+VecEE.y*Norm1.y), 
                                    fabs(VecEE.x*Norm2.x+VecEE.y*Norm2.y));
                            if (Elevation < 0.3){
                                VecE1.x = cos(Polygon[i].rot + ii*Angle*pi);
                                VecE1.y = sin(Polygon[i].rot + ii*Angle*pi);
                                VecE2.x = cos(Polygon[j].rot + jj*Angle*pi);
                                VecE2.y = sin(Polygon[j].rot + jj*Angle*pi);
                                Misalign = std::max(
                                        fabs(VecEE.x*VecE1.x + VecEE.y*VecE1.y),
                                        fabs(VecEE.x*VecE2.x + VecEE.y*VecE2.y));
                                if (Misalign <= EdgeLength){
                                    Overlap = 1 - Misalign/EdgeLength;
                                    NeighborF[i] += 1;
                                    NeighborF[j] += 1;
                                    LocF[i] += Overlap;
                                    LocF[j] += Overlap;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (int i = 0; i < np; i++){
        if (NeighborF[i] != 0){
            LocF[i] /= NeighborF[i];
            F += LocF[i];
        }
    }
    F /= np;
}
