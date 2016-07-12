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
void MeanSquaredDisplacement(std::string filename);

int main(int argc, char* argv[]) {
    if (argc == 1){
//	MeanSquaredDisplacement("traj/");
	std::cout << "Single thread quench simulation is starting.." << std::endl;
	testCppModel("traj/");
    } else if (argc == 2){
	std::cout << "Multi-thread quench simulation is starting.." << std::endl;
	testCppModelMT(boost::lexical_cast<int>(argv[1]));
//	std::cout << "Single thread reinforcement learning is starting.." << std::endl;
//    	testQLearning(argv[1]);
    } else if ( argc == 3) {
        int thread = boost::lexical_cast<int>(argv[2]);
        std::cout << thread << " thread reinforcement learning is starting.." << std::endl;
        testQLearningMT(argv[1],thread);
    } else {
	std::cout << "Wrong argument number" << std::endl;
    }
    return 0;
}

void testCppModel(std::string filename){
    int cycle(1),second(500);
    std::shared_ptr<BaseModel> model(new Model_QuadrupoleBD(filename));
    int iter;
    for (int i = 0; i < cycle; i++) {
        model->createInitialState();
        std::cout << "Cycle " << i << " is starting ..." << std::endl;
        iter = 0;
       while (iter < second) {
//        while (iter < second && !model->terminate()) {
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

void testQLearning(char* filename2){
    ReinforcementLearningParameter message2;
    QLearningSolverParameter message3;
    ReadProtoFromTextFile(filename2, &message2);
    message3 = message2.qlearningsolverparameter();
    std::shared_ptr<BaseModel> model(new Model_QuadrupoleBD("traj/control"));    

    int n_rows(5), n_cols(5);
    
    double dx1 = (double)1/n_rows;
    double dx2 = (double)1/n_cols;
    double minx1(0.0), minx2(0.0);
    
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
    
    int Resolution(5);
    int n_rows = Resolution;
    int n_cols = Resolution;
    
    double dx1 = (double) 1/n_rows;
    double dx2 = (double) 1/n_cols;
    double minx1(0.0), minx2(0.0);
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

void MeanSquaredDisplacement(std::string filename){
    std::shared_ptr<BaseModel> model(new Model_QuadrupoleBD(filename));
    std::vector<double> phi6start = {0.2, 0.4, 0.6, 0.8};
    int cycle = 0;
    while(!phi6start.empty()){
	double currstart = phi6start[0];
	model->createInitialState();
	std::cout << "cycle " << cycle << " target is Psi6 = " <<  currstart << std::endl;
	int second = 0;
	while (second < 100 && model->callpsi6() < currstart){
	    model -> run(3);
	    second++;
	    std::cout << "t = " << second << " psi6 = " << model->callpsi6() << std::endl;
	}
	if (model->callpsi6() >= currstart-0.005){
	    phi6start.erase(phi6start.begin());
	    std::cout << "psi6 = " << currstart << " achieved.." << std::endl;
	} else {
	    std::cout << "psi6 = " << currstart << " not achieved, start the next cycle..." << std::endl;
	}
	cycle++;
    }
}
