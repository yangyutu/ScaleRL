#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <memory>
#include "common_RL.h"
#include "Model_QuadrupoleBD.h"
#include "BaseModel.h"
#include "RLSolver_2DTable.h"
#include "RLSolver_2DTableMT.h"
#include <boost/lexical_cast.hpp>

using namespace ReinforcementLearning;

void testFortranModel();
void testCppModel(std::string filename);
void testCppModelMT(int n);
void readxyz(double *r, int np, const std::string filename);
void writexyz(double *r, int np, std::ostream &os);
void writeOp(double t, double psi6, double c6, double rg, double lambda, int opt, std::ostream &os);
void testQLearning(char* filename2);
void testQLearningMT(char* filename2, int t);


int main(int argc, char* argv[]) {
//    testCppModel("traj/");
//    testCppModelMT(boost::lexical_cast<int>(argv[1]));
//    testQLearning(argv[1]);
    if ( argc == 3) {
        std::cout << argc << " arguments" << std::endl;
        int thread = boost::lexical_cast<int>(argv[2]);
        std::cout << thread << " threads" << std::endl;
        testQLearningMT(argv[1],thread);
    }
    
    return 0;
}

void testCppModel(std::string filename){
    std::shared_ptr<BaseModel> model(new Model_QuadrupoleBD(filename));
    
    int iter;
    
    for (int i = 0; i < 1; i++) {
        model->createInitialState();
        std::cout << i << std::endl;
        iter = 0;
        while (iter < 10 && !model->terminate()) {
            model->run(3);
            iter++;
        }
    }    
}

void testCppModelMT(int nthreads){
    int num_threads = nthreads;
     std::thread *threads = new std::thread[num_threads];
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++){
        std::stringstream ss;
        ss << thread_idx;
        threads[thread_idx] = std::thread(testCppModel, "traj/" + ss.str());    
    }
    
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++){
        threads[thread_idx].join();    
    } 
    
     delete[] threads;
}

/*
void testFortranModel(){
    std::ofstream os1, os2;
    os1.open("bd_xyzC.txt");
    os2.open("opC.txt");
    double r[900];
    int n, opt, nstep;
    n = 300;
    nstep = 10000;
    opt = 3;
    double t, psi6, c6, rg, lambda;

    readxyz(r, n, "startmeshgrid1.txt");
    for (int i = 0; i < 10; i++) {
        run_fortran_(r, &n, &nstep, &psi6, &c6, &rg, &opt, &lambda);
        writexyz(r, n, os1);
        t = i * 0.1 * nstep / 1000.0;
        writeOp(t, psi6, c6, rg, lambda, opt, os2);
    }
}
*/
void readxyz(double *r, int np, const std::string filename) {
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
}

void writexyz(double *r, int np, std::ostream &os) {

    for (int i = 0; i < np; i++) {
        os << i << "\t";
        os << r[3 * i] << "\t";
        os << r[3 * i + 1] << "\t";
        os << r[3 * i + 2] << "\t";
        os << std::endl;
    }
}

void writeOp(double t, double psi6, double c6, double rg, double lambda, int opt, std::ostream &os) {
    os << t << "\t";
    os << psi6 << "\t";
    os << c6 << "\t";
    os << rg << "\t";
    os << opt << "\t";
    os << lambda << "\t";
    os << std::endl;
}

void testQLearning(char* filename2){
    ReinforcementLearningParameter message2;
    QLearningSolverParameter message3;
    ReadProtoFromTextFile(filename2, &message2);
    message3 = message2.qlearningsolverparameter();
    std::shared_ptr<BaseModel> model(new Model_QuadrupoleBD("traj/control"));
    
    int n_rows = 20;
    int n_cols = 20;
    
    double dx1 = 0.05;
    double dx2 = 0.05;
    double minx1 = 0.0;
    double minx2 = 0.0;
    
    RLSolver_2DTable rlSolver(model, 2, message3, n_rows, n_cols, dx1, dx2, minx1, minx2);
    rlSolver.loadQTable("QTableFinal");
//    rlSolver.getQTable().slice(2).ones();
    rlSolver.train();
}


void testQLearningMT(char* filename2, int thread){

    ReinforcementLearningParameter message2;
    QLearningSolverParameter message3;
    ReadProtoFromTextFile(filename2, &message2);
    message3 = message2.qlearningsolverparameter();
    
    int n_rows = 20;
    int n_cols = 20;
    
    double dx1 = 0.05;
    double dx2 = 0.05;
    double minx1 = 0;
    double minx2 = 0;
    
    int num_threads = thread;
    std::vector<std::shared_ptr<BaseModel>> models;
    for (int i = 0; i < num_threads; i++){
        std::stringstream ss;
        ss << i;
        models.push_back(std::shared_ptr<BaseModel>(new Model_QuadrupoleBD("traj/thread" + ss.str() + "control")));
    }
    
    
    RLSolver_2DTableMT rlSolver(models, 2, message3, n_rows, n_cols, dx1, dx2, minx1, minx2,num_threads);
    rlSolver.getQTable().slice(3).fill(5);
    rlSolver.train();
}