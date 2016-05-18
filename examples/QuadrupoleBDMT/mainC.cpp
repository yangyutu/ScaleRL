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
    if (argc == 1){
	testCppModel("traj/");
    } else if (argc == 2){
//    testCppModelMT(boost::lexical_cast<int>(argv[1]));
    	testQLearning(argv[1]);
    } else if ( argc == 3) {
        std::cout << argc << " arguments" << std::endl;
        int thread = boost::lexical_cast<int>(argv[2]);
        std::cout << thread << " threads" << std::endl;
        testQLearningMT(argv[1],thread);
    } else {
	std::cout << "argument number not recognized" << std::endl;
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
        std::cout << "t = " << iter << std::endl;
       while (iter < 10000) {
//        while (iter < 20 && !model->terminate()) {
            model->run(0);
            iter++;
            std::cout << "t = " << iter << std::endl;
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
    
    double dx1 = 1/20;
    double dx2 = 1/20;
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
    
    int Resolution = 5;
    int n_rows = 5;
    int n_cols = 5;
    
    double dx1 = 0.2;
    double dx2 = 0.2;
    double minx1 = 0.0;
    double minx2 = 0.0;
    int num_threads = thread;
    std::vector<std::shared_ptr<BaseModel>> models;
    for (int i = 0; i < num_threads; i++){
        std::stringstream ss;
        ss << i;
        models.push_back(std::shared_ptr<BaseModel>(new Model_QuadrupoleBD("traj/thread" + ss.str() + "control", Resolution)));
    }
    
    
    RLSolver_2DTableMT rlSolver(models, 2, message3, n_rows, n_cols, dx1, dx2, minx1, minx2,num_threads);
    rlSolver.getQTable().slice(3).fill(1);
//    rlSolver.loadQTable("./QTableFile/QTableFinal");
    rlSolver.train();
}
