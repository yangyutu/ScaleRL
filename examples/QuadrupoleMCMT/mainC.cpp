#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <memory>
#include "common_RL.h"
#include "Model_QuadrupoleMC.h"
#include "BaseModel.h"
#include "RLSolver_2DTable.h"
#include "RLSolver_2DTableMT.h"
#include <boost/lexical_cast.hpp>

using namespace ReinforcementLearning;

void Quench(std::string filename, int polygon, int thread_idx);
void RLMT(char* filename2, int thread, int polygon);
void RLST(char* filename2, int polygon);
void QuenchMT(int thread);

int main(int argc, char* argv[]) {
    if (argc == 1){
        std::cout << "Single Thread Quench Simulation is starting.." << std::endl;
        Quench("traj/",3,1);
    }else if ( argc == 2){
	int thread = boost::lexical_cast<int>(argv[1]);
	std::cout << thread << " Thread Quench Simulation is starting.." << std::endl;
        QuenchMT(thread);
    } else if ( argc == 3) {
	int polygon = boost::lexical_cast<int>(argv[2]);
	std::cout << "Single Thread Reinforcement Learning Simulation is starting.." << std::endl;
        RLST(argv[1],polygon);
    } else if ( argc == 4) {
        int thread = boost::lexical_cast<int>(argv[2]);
	int polygon = boost::lexical_cast<int>(argv[3]);
        std::cout << thread << " Thread Reinforcement Learning Simulation is starting.." << std::endl;
        RLMT(argv[1],thread,polygon);
    } else { 
        std::cout << "Argument number was not recognized" << std::endl;
    }
    return 0;
}

void Quench(std::string filename, int polygon, int thread_idx){
    int cycle(1), second(100);
    std::shared_ptr<BaseModel> model(new Model_QuadrupoleMC(filename,1,polygon));
    int iter;
    for (int i = 0; i < cycle; i++) {
        model->createInitialState();
        iter = 0;
//       while (iter < second) {
        while (iter < second && !model->terminate()) {
            model->run(3);
            iter++;
            std::cout << iter << std::endl;
        }
    std::cout << "Thread " << thread_idx+1 << " cycle " << i+1 << " completed" << std::endl;
    }
}

void QuenchMT(int nthreads){
    int num_threads = nthreads;
     std::thread *threads = new std::thread[num_threads];
    std::cout << "There are " << num_threads << " threads"<< std::endl;
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++){
        std::stringstream ss;
        ss << thread_idx;
        threads[thread_idx] = std::thread(Quench, "traj/" + ss.str(), 3, thread_idx);    
    }
    
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++){
        threads[thread_idx].join();   
    } 
    
     delete[] threads;
}

void RLST(char* filename2, int polygon){
    ReinforcementLearningParameter message2;
    QLearningSolverParameter message3;
    ReadProtoFromTextFile(filename2, &message2);
    message3 = message2.qlearningsolverparameter();

    int Resolution = 5;
    int n_rows = 5;
    int n_cols = 5;
    
    double dx1 = 1.0/5.0;
    double dx2 = 1.0/5.0;
    double minx1 = 0.0;
    double minx2 = 0.0;
    std::cout << "Starting a single thread RL simulation" << std::endl;
    std::shared_ptr<BaseModel> model(new Model_QuadrupoleMC("traj/control",Resolution,polygon));
    RLSolver_2DTable rlSolver(model, 2, message3, n_rows, n_cols, dx1, dx2, minx1, minx2);
//    rlSolver.loadQTable("QTableFinal");
    rlSolver.getQTable().slice(3).fill(1);
    rlSolver.train();
}

void RLMT(char* filename2, int thread, int polygon){

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
    std::cout << "Starting a " << thread << " thread RL simulation" << std::endl;
    std::vector<std::shared_ptr<BaseModel>> models;
    for (int i = 0; i < num_threads; i++){
        std::stringstream ss;
        ss << i;
        models.push_back(std::shared_ptr<BaseModel>(new Model_QuadrupoleMC("traj/thread" + ss.str() + "control", Resolution, polygon)));
    }
    
    
    RLSolver_2DTableMT rlSolver(models, 2, message3, n_rows, n_cols, dx1, dx2, minx1, minx2,num_threads);
    rlSolver.getQTable().slice(3).fill(1);
//    rlSolver.loadQTable("./QTableFile/QTableFinal");
    rlSolver.train();
}
