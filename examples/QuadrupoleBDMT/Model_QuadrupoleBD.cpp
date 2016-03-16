#include "Model_QuadrupoleBD.h"
using namespace ReinforcementLearning;
// this model is from paper Lease-squares policy iteration

Model_QuadrupoleBD::Model_QuadrupoleBD(std::string filetag0) {
    filetag = filetag0;
    //	we have three low dimensional states psi6, c6, rg
    // nstep 10000 correspond to 1s, every run will run 1s
    nstep = 10000;
    stateDim = 3;
    currState.resize(stateDim);
    prevState.resize(stateDim);
    dt = 1; //s
    numActions = 4;
    trajOutputInterval = 1;
    fileCounter = 0;
    rand_normal = std::make_shared<std::normal_distribution<double>>(0.0, 1.0);
    rand_int = std::make_shared<std::uniform_int_distribution<>>(0, 1000);
    rgbin = 25;
    rbin = 50;
    a = 1435.0;
    kb = 1.380658e-23;
 


    for (int i = 0; i < np; i++) {
        nxyz[i][0] = 3 * i;
        nxyz[i][1] = 3 * i + 1;
        nxyz[i][2] = 3 * i + 2;
    }
}

void Model_QuadrupoleBD::run(int action) {
    this->opt = action;
    if (this->timeCounter == 0 || ((this->timeCounter + 1) % trajOutputInterval == 0)) {
        this->outputTrajectory(this->trajOs);
        this->outputOrderParameter(this->opOs);
    }

    int rand = (*rand_int)(rand_generator);
    rand = 1;
    // for BD dynamics, every run will simply run 10000 steps, correspond to 1s
    //run_fortran_(r, &np, &nstep, &psi6, &c6, &rg, &opt, &lambda, &rand, dss, dssCount);
    this->runHelper(nstep,opt);
    
    c6 = c6 / 5.6;
    this->prevState = this->currState;
    this->currState[0] = psi6;
    this->currState[1] = c6;
    this->timeCounter++;
}

void Model_QuadrupoleBD::createInitialState() {
    if (fileCounter < 15) {this->readxyz("./StartMeshgridFolder/startmeshgrid1.txt");}
    else {
    std::stringstream FileStr;
    FileStr << this->fileCounter;
    this->readxyz("./StartMeshgridFolder/startmeshgrid" + FileStr.str() + ".txt");
    }
    this->readDiffusivity("2dtabledsslam9.txt");
    std::stringstream ss;
    std::cout << "model initialize at round " << fileCounter << std::endl;
    ss << this->fileCounter++;
    if (trajOs.is_open()) trajOs.close();
    if (opOs.is_open()) opOs.close();

    this->trajOs.open(filetag + "xyz_" + ss.str() + ".dat");
    this->opOs.open(filetag + "op" + ss.str() + ".dat");
    this->timeCounter = 0;

    int temp = 0;
    opt = 0;
    int rand = (*rand_int)(rand_generator);
    rand = 1;
    this->runHelper(temp,opt);
    //run_fortran_(r, &np, &temp, &psi6, &c6, &rg, &opt, &lambda, &rand, dss, dssCount);
    c6 = c6 / 5.6;
    this->currState[0] = psi6;
    this->currState[1] = c6;

}

void Model_QuadrupoleBD::outputTrajectory(std::ostream& os) {

    for (int i = 0; i < np; i++) {
        os << i << "\t";
        os << r[3 * i]/a << "\t";
        os << r[3 * i + 1]/a << "\t";
        os << r[3 * i + 2]/a<< "\t";
        os << std::endl;
    }
}

void Model_QuadrupoleBD::outputOrderParameter(std::ostream& os) {
    os << this->timeCounter * dt << "\t";
    os << psi6 << "\t";
    os << c6 << "\t";
    os << rg << "\t";
    os << opt << "\t";
    os << lambda << "\t";
    os << std::endl;
}

double Model_QuadrupoleBD::getRewards() {
    if (this->terminate()) {
        this->reward = 0;
        return reward;
    } else {
        this->reward = -(1 - currState[0]);
        return reward;
    }
}

bool Model_QuadrupoleBD::terminate() {
    // the model will stop if \psi6 > 0.95, which is perfect crystal
    if (currState[0] > 0.90) {
        return true;
    }
    return false;
}

void Model_QuadrupoleBD::readxyz(const std::string filename) {
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
        linestream >> r[3 * i + 2];
    }
    for (int i = 0; i < np * 3; i++) {
        r[i] *= a;
    }

    is.close();
}

void Model_QuadrupoleBD::readDiffusivity(const std::string filename) {
    std::ifstream is;
    is.open(filename.c_str());
    std::string line;
    double dum;
    for (int i = 0; i < rbin * rgbin; i++) {
        getline(is, line);
        std::stringstream linestream(line);
        linestream >> dum;
        linestream >> dum;
        linestream >> dum;
        linestream >> dum;
        linestream >> *(&(dssarray[0][0])+i);
        linestream >> *(&(dsscount[0][0])+i);
    }
    is.close();
}

void Model_QuadrupoleBD::runHelper(int nstep, int controlOpt) {

    tempr = 20.0;
    double dt = 0.1;
    fac1 = 5.9582e7;
    fac2 = 40.5622;
    rcut = 5.0 * a;
    re = 5.0 * a;
    kappa = 143.5;
    pfpp = 2.2975 * a;
    fcm = -0.4667;
    DG = 71.428 * a;
    rmin = 3780;
    rgdsmin = 22250;
    delrgdsmin = -250;
    distmin = 0;
    deldist = 1400;
    dssmin = 0.15;
    dssmax = 0.5;
    dpf = 1;
    ecorrectflag = 1;


    if (controlOpt == 0) {
        lambda = 0.219;
    } else if (controlOpt == 1) {
        lambda = 0.8744;
    } else if (controlOpt == 2) {
        lambda = 1.9674;
    } else {
        lambda = 19.73;
    }


    fac1 = fac1 / a;
    fac2 = fac2 * sqrt((273 + tempr) / a);
    fac2 = fac2 / sqrt(dt);

    int step = 0;
    this->calDss();
    while (step < nstep) {

        
        for (int j = 0; j < np; j++) {
            for (int k = 0; k < 3; k++) {
                D[nxyz[j][k]] = dsscalcu[j];
                double randTemp = (*rand_normal)(rand_generator);
                randisp[nxyz[j][k]] =  randTemp * sqrt(1.0 / D[nxyz[j][k]]);
            }
        }

        forces();
        double u;
        for (int j = 0; j < np3; j++) {
            u = D[j] * (F[j] * fac1 + randisp[j] * fac2);
            r[j] += u * dt;
        }
        step++;
    }
    this->calOp();
}

void Model_QuadrupoleBD::forces() {
    double RX, RY, EMAGI, EMAGJ, dE2x, dE2y, Fdepx, Fdepy;
    double STEP = 1e-3;
    double rij[3];
    double Fpp, Fhw, felx, fely, felz, Exi, Exj, Eyi,Eyj,Ezi,Ezj;
    double felxnew,felxnew2,felynew,felynew2;
    Fhw = 0.417;
    double Fo = 1e18*0.75*lambda*kb*(273+tempr)/a;
    for (int i = 0 ; i < np3; i++) {
        F[i] = 0.0;
    }
    for (int i = 0; i < np-1; i++){
        Exi = -4.0*r[nxyz[i][0]]/DG;
        Eyi = 4.0*r[nxyz[i][1]]/DG;
        Ezi = 0;
        
        
        for (int j = i+1; j < np; j++){
            rij[0] = r[nxyz[j][0]] - r[nxyz[i][0]];
            rij[1] = r[nxyz[j][1]] - r[nxyz[i][1]];
            
            double rijsep = sqrt(rij[0]*rij[0]+rij[1]*rij[1]);
            if (rijsep < 2*a){
                Fpp = Fhw;
                felx = 0.0;
                fely = 0.0;
                felz = 0.0;
            } else if( rijsep < rcut) {
                Fpp = 1e18*kb*(tempr+273)*kappa*pfpp*exp(-kappa*(rijsep - 2.0*a)/a)/a;
            } else {
                Fpp = 0;
            }
            
            if (rijsep > 2*a && rijsep < re){
                Exj = -4.0*r[nxyz[j][0]]/DG;
                Eyj = 4.0*r[nxyz[j][1]]/DG;
                Ezj = 0;
             
                double F1 = Exi*Exj + Eyi*Eyj + Ezi*Ezj;
                double F2 = rij[0]*Exi/rijsep + rij[1]*Eyi/rijsep + rij[2]*Ezi/rijsep;
                double F3 = rij[0]*Exj/rijsep + rij[1]*Eyj/rijsep + rij[2]*Ezj/rijsep;
                
                felxnew=Fo*pow(2*a/rijsep,4)*(F1*rij[0]/rijsep + Exi*F3 + 
                        Exj*F2 - 5*F2*F3*rij[0]/rijsep+
                        rijsep*16*r[nxyz[j][0]]/pow(DG,2)/3.0-F3*4*(-rij[0])/DG);
                felynew=Fo*pow(2*a/rijsep,4)*(F1*rij[1]/rijsep + Eyi*F3 + 
                        Eyj*F2 - 5*F2*F3*rij[1]/rijsep+
                        rijsep*16*r[nxyz[j][1]]/pow(DG,2)/3.0-F3*4*(rij[1])/DG);
                
                
                     
                felxnew2=Fo*pow(2*a/rijsep,4)*(F1*rij[0]/rijsep +
                        Exi*F3 + Exj*F2 - 5*F2*F3*rij[0]/rijsep-
                        rijsep*16*r[nxyz[i][0]]/pow(DG,2)/3.0+F2*4*(-rij[0])/DG);
                felynew2=Fo*pow(2*a/rijsep,4)*(F1*rij[1]/rijsep +
                        Eyi*F3 + Eyj*F2 - 5*F2*F3*rij[1]/rijsep-
                        rijsep*16*r[nxyz[i][1]]/pow(DG,2)/3.0+F2*4*(rij[1])/DG);
                       
            } else {

                felx = 0.0;
		fely = 0.0;
		felz = 0.0;
		felxnew=0.0;
		felynew=0.0;
		felxnew2=0.0;
		felynew2=0.0;
            
            }
            
            
            
            F[nxyz[i][0]] = F[nxyz[i][0]] - felxnew - Fpp*rij[0]/rijsep;
            F[nxyz[i][1]] = F[nxyz[i][1]] - felynew - Fpp*rij[1]/rijsep;
            
            F[nxyz[j][0]] = F[nxyz[j][0]] + felxnew2 + Fpp*rij[0]/rijsep;
            F[nxyz[j][1]] = F[nxyz[j][1]] + felynew2 + Fpp*rij[1]/rijsep;
            
        }
    }
       
    for (int i = 0; i < np; i++) {
        RX = r[nxyz[i][0]];
        RY = r[nxyz[i][1]];
        EMAGI = EMAG(RX, RY);
        
        RX = r[nxyz[i][0]] + STEP;
        RY = r[nxyz[i][1]];
        EMAGJ = EMAG(RX, RY);
        
        dE2x = (EMAGJ * EMAGJ - EMAGI * EMAGI) / STEP;
        RX = r[nxyz[i][0]];
        RY = r[nxyz[i][1]] + STEP;

        EMAGJ = EMAG(RX, RY);
        dE2y = (EMAGJ * EMAGJ - EMAGI * EMAGI) / STEP;
        Fdepx = (2 * 1e18 * kb * (tempr + 273) * lambda / fcm) * dE2x;
        Fdepy = (2 * 1e18 * kb * (tempr + 273) * lambda / fcm) * dE2y;
        F[nxyz[i][0]] += Fdepx;
        F[nxyz[i][1]] += Fdepy;

    }
}

double Model_QuadrupoleBD::EMAG(double RX, double RY) {

    double correctfactor;

    double RT = sqrt(RX * RX + RY * RY);
    double result = 4 * RT / DG;
    if (ecorrectflag == 1) {
        correctfactor = 2.081e-7 * pow(RT / 1000.0, 4) - 1.539e-9 * pow(RT / 1000.0, 3) + 8.341e-5 * pow(RT / 1000, 2) + 1.961e-5 * (RT / 1000) + 1.028;
        result *= correctfactor;
    }
    return result;

}

void Model_QuadrupoleBD::calOp() {

    int nb[np], con[np];
    double rx[np], ry[np];
    double rxij, ryij, theta, psir[np], psii[np], numer, denom, testv, ctestv;
    double rgmean, xmean, ymean, accumpsi6r, accumpsi6i;
    ctestv = 0.32;
    for (int i = 0; i < np; i++) {
        rx[i] = r[nxyz[i][0]];
        ry[i] = r[nxyz[i][1]];
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

void Model_QuadrupoleBD::calDss() {

    int rgbinindex = (int) ((rg - rgdsmin) / delrgdsmin) + 1;

    if (rgbinindex <= 0) {
        rgbinindex = 1;
    }

    double calcudss = 0.5 * (dssmax + dssmin);
    double xmean = 0.0;
    double ymean = 0.0;

    for (int i = 0; i < np; i++) {
        xmean = xmean + r[nxyz[i][0]];
        ymean = ymean + r[nxyz[i][1]];
    }
    xmean = xmean / np;
    ymean = ymean / np;

    for (int i = 0; i < np; i++) {
        double disttemp = pow(r[nxyz[i][0]] - xmean, 2) + pow(r[nxyz[i][1]] - ymean, 2);
        disttemp = sqrt(disttemp);
        int distbinindex = (int) ((disttemp - distmin) / deldist) + 1;
        if (rgbinindex >= 1 && rgbinindex <= rgdssbin) {
            if (distbinindex >= 1 && distbinindex <= distdssbin) {

                if (dsscount[rgbinindex-1][distbinindex-1] >= 1) {
                    dsscalcu[i] = dssarray[rgbinindex - 1][distbinindex - 1];
                } else {
                    dsscalcu[i] = dssmax;
                }
            } else {

                dsscalcu[i] = dssmax;
            } 
        } else if (rgbinindex >= rgdssbin) {
            dsscalcu[i] = dssmin;
        } 

        }
}
