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
void writeCountMap(arma::cube CountMap);
void RLMT(char* filename2, int thread, int polygon);

int main(int argc, char* argv[]) {
    if ( argc == 2){
	int polygon = boost::lexical_cast<int>(argv[1]);
	Quench("traj/",polygon);
    } else if ( argc == 4) {
        std::cout << argc << " arguments" << std::endl;
        int thread = boost::lexical_cast<int>(argv[2]);
	int polygon = boost::lexical_cast<int>(argv[3]);
        std::cout << thread << " threads" << std::endl;
        RLMT(argv[1],thread,polygon);
    } else { std::cout << "Argument number was not recognized" << std::endl;}
    return 0;
}

void Quench(std::string filename, int polygon){
    std::shared_ptr<BaseModel> model(new Model_QuadrupoleMC(filename,1,polygon));
    int iter;
    for (int i = 0; i < 1; i++) {
        model->createInitialState();
        iter = 0;
       while (iter < 20) {
//        while (iter < 300 && !model->terminate()) {
            model->run(3);
	    std::cout << "t = " << iter << std::endl;
            iter++;
        }
    }
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

void writeCountMap(arma::cube CountMap) {
    arma::mat temp1 = CountMap.slice(1);
    temp1.save("Option1.dat", arma::raw_ascii);
    arma::mat temp2 = CountMap.slice(2);
    temp2.save("Option2.dat", arma::raw_ascii);
    arma::mat temp3 = CountMap.slice(3);
    temp3.save("Option3.dat", arma::raw_ascii);
    arma::mat temp4 = CountMap.slice(4);
    temp4.save("Option4.dat", arma::raw_ascii);
}
