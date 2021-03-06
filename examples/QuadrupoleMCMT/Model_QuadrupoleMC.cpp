#include "Model_QuadrupoleMC.h"

using namespace ReinforcementLearning;

Model_QuadrupoleMC::Model_QuadrupoleMC(std::string filetag0,int Resolution, int polygon0) {
//    File information
    filetag = filetag0;     // Path of output files
    polygonnum = polygon0;  // Number of polygon edges

    //    Physical constants
    DiffTrans = 6.362e-5;   // Diffusivity of translational motion
    DiffRot = 4.772e-5;     // Diffusivity of rotational motion
    Angle = 2.0/polygonnum; // Degree of symmetry

    // 1 is the length from center to corner
    a = sin(0.5*pi-pi/polygonnum);      // the length from center to edge
    EdgeLength = a*tan(pi/polygonnum);  // the length of half of edge
    rmin = 2.2;             // Distance criterion for Psi6
    rmin2 = 3*a;            // Distance criterion for F
    rmin3 = 2.2*a;          // Distance criterion for Chi
    ctestv = 0.32;          // Distance criterion for C6
    epsilon = 8.854e-12*80; // Dielectric permittivity of water [SI]
    kT = 1.38e-23*293;      // Boltzmann Constant * temperature [SI]
    e = 1.6e-19;            // Elementary charge [C]
    z = 50e-3;              // Surface charge potential [V]
    fcm = -0.4667;
    dg = 100/1.5;
    OSM = 0.0149;

    //    Simulation parameters
    nstep = 10000;          // Number of steps recurse in each time intervial (1s)
    dt = 1000.0/nstep;      // Time of each step recursed (1ms)
    stateDim = 5;           // State of polygon configuration in Ps6, Rg, F
    currState.resize(stateDim);
    prevState.resize(stateDim);
    numActions = 4;         // Number of voltages used
    fileCounter = 0;        // Cycle number in each thread
    n_rows = Resolution;
    n_cols = Resolution;
    dx1 = 1.0/Resolution;
    dx2 = 1.0/Resolution;
    rand_normal = std::make_shared<std::normal_distribution<double>>(0.0, 1.0);
    rand_uniform = std::uniform_real_distribution<double> (0, 1);
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
    coord temp;
    for (int i = 0; i < np; i++){ // Coordinate of edges and particle zone
        Polygon[i].edge.clear();
        for (int j = 0; j < polygonnum; j++){
            temp.x = Polygon[i].center.x + a*cos(Polygon[i].rot - pi/2 + j*Angle*pi);
            temp.y = Polygon[i].center.y + a*sin(Polygon[i].rot - pi/2 + j*Angle*pi);
            Polygon[i].edge.push_back(temp);
        }
    }
    for (int i = 0; i < np; i++){
        for (int j = i + 1; j < np; j++){
            if(Distance(Polygon[i],Polygon[j]) <= 5*a){
                Polygon[i].Neighbor.push_back(j);
                Polygon[j].Neighbor.push_back(i);
            }
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
    this->runCore(0,3);
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
//        linestream >> dum;
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
	os << psi6 << "\t";
	os << rg << "\t";
        os << Polygon[i].F << "\t";
        os << Polygon[i].C << "\t";
        os << std::endl;
    }
}

void Model_QuadrupoleMC::outputOrderParameter(std::ostream& os) {
    os << this->timeCounter << "\t";
    os << psi6 << "\t";
    os << rg << "\t";
    os << F << "\t";
    os << C << "\t";
    os << opt << "\t";
    os << lambda << "\t";
    os << std::endl;
}

double Model_QuadrupoleMC::getRewards() {
    if (this->terminate()) {
        this->reward = 0;
        return reward;
    } else {
        this->reward = -(1 - (currState[0]+currState[1])/2.0);
        return reward;
    }
}

bool Model_QuadrupoleMC::terminate() {
    if (currState[0] > 0.90 && currState[1] > 0.95) {return true;}
    return false;
}

void Model_QuadrupoleMC::run(int action) {
/* Takes in the control option and run the dynamics for 1s*/
    this->opt = action;
    this->prevState = this->currState;
    this->runCore(nstep,opt);
    this->timeCounter++;
}

void Model_QuadrupoleMC::runCore(int nstep, int controlOpt) {
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
    // Run 1s of simulation in nstep steps; update the neighbor every 100 loops
    for (int step = 0; step < nstep; step++) {
        this->MonteCarlo();
        if ((step+1)%100 == 0){
            UpdateNeighbor();
        }
    }
    // Output order parameters and trajectories
    this->calPsi();
    this->calRg();
    this->calF();
    this->calC();
    this->currState[0] = F;
    this->currState[1] = psi6;
    this->currState[2] = rg;
    this->currState[3] = C;
    this->outputTrajectory(this->trajOs);
    this->outputOrderParameter(this->opOs);
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
 *      movement is rejected following a Boltzmann distribution. 
 */
    
    double Driftx, Drifty, RandDriftx, RandDrifty, RandRot;
    for (int i = 0; i < np; i++){
// The velocity of each particle (translational and rotational)
        particle TempPolygon = Polygon[i];
        Driftx = -Polygon[i].center.x*lambda*dt*DiffTrans;
        Drifty = -Polygon[i].center.y*lambda*dt*DiffTrans;
        RandDriftx = (*rand_normal)(rand_generator);
        RandDrifty = (*rand_normal)(rand_generator);
        RandRot = (*rand_normal)(rand_generator);

// Update new structure in Polygon
        Polygon[i].center.x += Driftx + RandDriftx*sqrt(DiffTrans*2.0*dt);
        Polygon[i].center.y += Drifty + RandDrifty*sqrt(DiffTrans*2.0*dt);
        Polygon[i].rot += RandRot*sqrt(DiffRot*2.0*dt);
        for (int edge = 0; edge < polygonnum; edge++){
            Polygon[i].edge.at(edge).x = Polygon[i].center.x + a*cos(Polygon[i].rot - pi/2 + edge*Angle*pi);
            Polygon[i].edge.at(edge).y = Polygon[i].center.y + a*sin(Polygon[i].rot - pi/2 + edge*Angle*pi);
        }
        
        // Check overlap and energy change
        bool OverLap = false;   // true = overlap occurs 
        bool Energy_pref = true;       // true = energy prefers movement    
        // Check particles that are within 4 blocks from particle i
        for (std::vector<int>::iterator j = Polygon[i].Neighbor.begin(); j != Polygon[i].Neighbor.end(); ++j){
            double tempdist = Distance(Polygon[i],Polygon[*j]);   
            if (tempdist < 2*a){
                OverLap = true;
            }
            if (this->CheckOverlap(Polygon[i],Polygon[*j]) == true){
                OverLap = true;
            }                
        }
// If overlap, the movement is retrieved; otherwise, the movement is accepted at
// the probability of Boltzmann distribution
        if ( OverLap == false){
            Energy_pref = CheckEnergy(Polygon[i], TempPolygon);
            if (Energy_pref == false){ // retrieve move if energy unfavored
                Polygon[i] = TempPolygon;
            }
        } else {
            Polygon[i] = TempPolygon;
        }
    }
}

inline bool Model_QuadrupoleMC::CheckOverlap(particle i, particle j){
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
                    Det1 = -Diffx*sin(j.rot + jj*Angle*pi) + 
                            Diffy*cos(j.rot + jj*Angle*pi);
                    Det2 = Diffy*cos(i.rot + ii*Angle*pi) - 
                            Diffx*sin(i.rot + ii*Angle*pi); 
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

inline bool Model_QuadrupoleMC::CheckEnergy(particle Pnew, particle Pold){
    double EnergyDiff, Boltzmann, RanGen;
    double U_pf_old, U_pf_new;
    double U_dp_old, U_dp_new;
    U_pf_old = 2*kT*lambda/fcm*4*sqrt(pow(Pold.center.x,2.0) + pow(Pold.center.y,2.0))/dg;
    U_pf_new = 2*kT*lambda/fcm*4*sqrt(pow(Pnew.center.x,2.0) + pow(Pnew.center.y,2.0))/dg;
    U_dp_old = 0;
    U_dp_new = 0;
    for (std::vector<int>::iterator j = Pold.Neighbor.begin(); j != Pold.Neighbor.end(); j++){
        double r_old = Distance(Pold, Polygon[*j]);
        double r_new = Distance(Pnew, Polygon[*j]);
        if (r_old <= 1.2*a){
            U_dp_old += -(4/3)*OSM*pi*1.2*a*pow(kT/e,3)*(1-0.75*r_old/(1.2*a)+0.0625*pow(r_old,3)/pow(1.2*a,3));
        }
        if (r_new <= 1.2*a){
            U_dp_new += -(4/3)*OSM*pi*1.2*a*pow(kT/e,3)*(1-0.75*r_new/(1.2*a)+0.0625*pow(r_new,3)/pow(1.2*a,3));
        }
    }
    EnergyDiff = (U_pf_new + U_dp_new) - (U_pf_old + U_dp_old);
    if (EnergyDiff > 0)
    {
        Boltzmann = exp(-EnergyDiff/kT);
        RanGen = rand_uniform(rand_generator);
        if (RanGen <= Boltzmann)
        { 
            return true;
        } 
        else 
        {
            return false;
        }
    } else 
    { 
        return true;
    }
}

void Model_QuadrupoleMC::calPsi() {

    int nb[np];
    double psir[np], psii[np], scale;
    double accumpsi6r, accumpsi6i;
    
    switch (polygonnum)
    {
        case 3:
            scale = 6.0;
            break;
        case 4:
            scale = 4.0;
            break;
        case 6:
            scale = 6.0;
            break;
        default:
            scale = 1.0;
            break;
    }
// Initialization
    psi6 = 0;
    accumpsi6r = 0;
    accumpsi6i = 0;
    for (int i = 0; i < np; i++) {
        nb[i] = 0;
        psir[i] = 0.0;
        psii[i] = 0.0;
    }
// calculate local psi6 in complex form
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < np; j++) {
            if (i != j) {
                double RP = sqrt(
                pow((Polygon[i].center.x-Polygon[j].center.x),2.0) + 
                pow((Polygon[i].center.y-Polygon[j].center.y),2.0));
                if (RP <= rmin) {
                    nb[i] += 1;
                    psir[i] += cos(scale * Polygon[j].rot);
                    psii[i] += sin(scale * Polygon[j].rot);
                }
            }        
        }
    }
    for (int i = 0; i < np; i++){
        if (nb[i] > 0) {
            psir[i] /=  nb[i];
            psii[i] /=  nb[i];
        }
    }
// Calculate global psi
    for (int i = 0; i < np; i++) {
        accumpsi6r += psir[i];
        accumpsi6i += psii[i];
    }
    accumpsi6r /= np;
    accumpsi6i /= np;
    psi6 = sqrt(accumpsi6r * accumpsi6r + accumpsi6i * accumpsi6i);
}

void Model_QuadrupoleMC::calRg() {
    double rgmean, xmean, ymean;
// Initialization
    xmean = 0;
    ymean = 0;
    rgmean = 0;
// Calculate mean coordinate
    for (int i = 0; i < np; i++) {
        xmean += Polygon[i].center.x;
        ymean += Polygon[i].center.y;
    }
    xmean /= np;
    ymean /= np;
// Calculate local rg
    for (int i = 0; i < np; i++) {
         Polygon[i].Rg = sqrt(pow(
                (Polygon[i].center.x - xmean),2) + 
                pow((Polygon[i].center.y - ymean),2));
         rgmean += (Polygon[i].Rg*Polygon[i].Rg)/np;
        
    }
// Calculate global rg
    rg = sqrt(rgmean);
}

void Model_QuadrupoleMC::calF() {
    /* New Order Parameter F.
     * 1.   The center-center distance must be shorter than 2.5a
     * 2.   Particle 2 must "face" an edge of particle 1 
     * 3.   Particle 2 must has an edge close to the edge of particle 1
     */
    coord VecCC, Norm1, VecEE, Norm2, VecE1, VecE2;
    double DistCC,angle, Elevation, Misalign, Overlap, Interception;
    int NeighborF[np];
    
    switch (polygonnum)
    {
        case 3:
             Interception = 1/(sqrt(7)/2);
	     break;
        case 4:
	    Interception = sqrt(2.0)/(sqrt(5.0)/sqrt(2.0));
	    break;
        case 6:
            Interception = sqrt(3.0)/(sqrt(13)/2);
	     break;
        default:
            Interception = sqrt(3.0)/2.0;
	     break;
    }
// Initialization
    F = 0;
    for (int i = 0; i < np; i++){
        Polygon[i].F = 0.0;
        NeighborF[i] = 0;        
    }
// Calculate local F
    for (int i = 0; i < np; i++){
        for (int j = i+1; j < np; j++){
            // Center-center distance of two particles, if close enough, calculate overlap
            DistCC = sqrt(pow((Polygon[i].center.x-Polygon[j].center.x),2.0) 
                    + pow((Polygon[i].center.y-Polygon[j].center.y),2.0));
            if (DistCC < rmin2){
                // Vector connecting two particle centers
                VecCC.x = (Polygon[j].center.x-Polygon[i].center.x);
                VecCC.y = (Polygon[j].center.y-Polygon[i].center.y);
                for (int ii = 0; ii < polygonnum; ii++){
                    // normal vector of each edge of particle 1
                    Norm1.x = (Polygon[i].edge.at(ii).x - Polygon[i].center.x)/a;
                    Norm1.y = (Polygon[i].edge.at(ii).y - Polygon[i].center.y)/a;
                    // Angle between Norm1 and VecCC
                    angle = (Norm1.x*VecCC.x + Norm1.y*VecCC.y)
                            /sqrt(VecCC.x*VecCC.x + VecCC.y*VecCC.y);
                    if (angle >= Interception){
                        for (int jj = 0; jj < polygonnum; jj++){
                            // Edge - edge vector between all edges
                            VecEE.x = Polygon[j].edge.at(jj).x - Polygon[i].edge.at(ii).x;
                            VecEE.y = Polygon[j].edge.at(jj).y - Polygon[i].edge.at(ii).y;
                            // normal vector of each edge of particle 2
                            Norm2.x = (Polygon[j].edge.at(jj).x - Polygon[j].center.x)/a;
                            Norm2.y = (Polygon[j].edge.at(jj).y - Polygon[j].center.y)/a;
                            // Vertical distance between two edges
                            Elevation = std::max(
                                    fabs(VecEE.x*Norm1.x+VecEE.y*Norm1.y), 
                                    fabs(VecEE.x*Norm2.x+VecEE.y*Norm2.y));
                            if (Elevation < a){
                                // Parallel vectors of each edge of each particle                 
                                VecE1.x = cos(Polygon[i].rot + ii*Angle*pi);
                                VecE1.y = sin(Polygon[i].rot + ii*Angle*pi);
                                VecE2.x = cos(Polygon[j].rot + jj*Angle*pi);
                                VecE2.y = sin(Polygon[j].rot + jj*Angle*pi);
                                // distance between edge center on edge directions
                                Misalign = std::max(
                                        fabs(VecEE.x*VecE1.x + VecEE.y*VecE1.y),
                                        fabs(VecEE.x*VecE2.x + VecEE.y*VecE2.y));
                                if (Misalign <= 1.2*EdgeLength){
                                    Overlap = 1 - Misalign/EdgeLength;
                                    if (Overlap <= 0){
                                        Overlap = 0;
                                    }
                                    NeighborF[i] += 1;
                                    NeighborF[j] += 1;
                                    Polygon[i].F += Overlap;
                                    Polygon[j].F += Overlap;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
// Calculate global F
    for (int i = 0; i < np; i++){
        if (NeighborF[i] != 0){
            Polygon[i].F /= NeighborF[i];
            F += Polygon[i].F;
        }
    }
    F /= np;
}

void Model_QuadrupoleMC::calC() {

    int neighbor[np];
    double psir[np], psii[np], scale;
    double numerator, denominator,testv;
    
    switch (polygonnum)
    {
        case 3:
            scale = 6.0;
            break;
        case 4:
            scale = 4.0;
            break;
        case 6:
            scale = 6.0;
            break;
        default:
            scale = 1.0;
            break;
    }
// Initialize OP
    C = 0.0;
    for (int i = 0; i < np; i++){
        neighbor[i] = 0;
        psir[i] = 0.0;
        psii[i] = 0.0;
        Polygon[i].C = 0.0;
    }
// Calculate local psi6 in complex form
    for (int i = 0; i < np; i++) {

        for (int j = 0; j < np; j++) {
            if (i != j) {
                double RP = sqrt(
                pow((Polygon[i].center.x-Polygon[j].center.x),2.0) + 
                pow((Polygon[i].center.y-Polygon[j].center.y),2.0));
                if (RP <= rmin) {
                    neighbor[i] += 1;
                    psir[i] += cos(scale * Polygon[j].rot);
                    psii[i] += sin(scale * Polygon[j].rot);
                }
            }        
        } 
        if (neighbor[i] > 0) {
            psir[i] /=  neighbor[i];
            psii[i] /=  neighbor[i];
        }
    }
// Calculate local C6
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < np; j++) {
            double RP = sqrt(
                pow((Polygon[i].center.x-Polygon[j].center.x),2.0) + 
                pow((Polygon[i].center.y-Polygon[j].center.y),2.0));
            if ((i != j)&& (RP <= rmin)) {
                numerator = psir[i] * psir[j] + psii[i] * psii[j];
                double temp = psii[i] * psir[j] - psii[j] * psir[i];
                denominator = sqrt(numerator * numerator + temp*temp);
                testv = numerator / denominator;
                if (testv >= ctestv) {
                    Polygon[i].C += 1;
                }
            }
        }
    }
// Average local C to get global C6
    for (int i = 0; i < np; i++){
        Polygon[i].C /= 10;
        C += Polygon[i].C;
    }
    C /= np;
}

void Model_QuadrupoleMC::UpdateNeighbor(){
    for (int i = 0; i < np; i++){
        Polygon[i].Neighbor.clear();
    }
    for (int i = 0; i < np; i++){
        for (int j = i + 1; j < np; j++){
            if(Distance(Polygon[i],Polygon[j]) < 5*a){
                Polygon[i].Neighbor.push_back(j);
                Polygon[j].Neighbor.push_back(i);
            }
        }
    }
}

double Model_QuadrupoleMC::Distance(particle i, particle j){
    double dx = i.center.x - j.center.x;
    double dy = i.center.y - j.center.y;
    return sqrt(dx*dx + dy*dy);
}