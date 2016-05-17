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

void Quench(std::string filename, int polygon);
void RLMT(char* filename2, int thread, int polygon);
void RLST(char* filename2, int polygon);

int main(int argc, char* argv[]) {
    if ( argc == 2){
	int polygon = boost::lexical_cast<int>(argv[1]);
	Quench("traj/",polygon);
    } else if ( argc == 3) {
	int polygon = boost::lexical_cast<int>(argv[2]);
	RLST(argv[1],polygon);
    } else if ( argc == 4) {
        int thread = boost::lexical_cast<int>(argv[2]);
	int polygon = boost::lexical_cast<int>(argv[3]);
        RLMT(argv[1],thread,polygon);
    } else { 
        std::cout << "Argument number was not recognized" << std::endl;
    }
    return 0;
}

void Quench(std::string filename, int polygon){
    std::cout << "Starting a quench sim" << std::endl;
    std::shared_ptr<BaseModel> model(new Model_QuadrupoleMC(filename,1,polygon));
    int iter;
    for (int i = 0; i < 1; i++) {
        model->createInitialState();
        iter = 0;
        std::cout << "t = 0" << std::endl;
       while (iter < 100) {
//        while (iter < 300 && !model->terminate()) {
            model->run(3);
	    std::cout << "t = " << iter+1 << std::endl;
            iter++;
        }
    }
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
    std::cout << "Starting a single thread RL sim" << std::endl;
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
    std::cout << "Starting a " << thread << " thread RL sim" << std::endl;
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
