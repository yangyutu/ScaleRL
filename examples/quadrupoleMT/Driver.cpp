#include <iostream>
#include "boost/multi_array.hpp"
#include "common_RL.h"
#include "Model_QuadrupoleLD.h"
#include "RLSolver_2DTableMT.h"
using namespace ReinforcementLearning;

void testModel();
void testQLearningMT(char* filename);
void testQLearning(char* filename);
int main(int argc, char* argv[]){
//    testModel();
//    testQLearning(argv[1]);
    testQLearningMT(argv[1]);
    return 1;
}


void testModel() {
    State state;
    std::shared_ptr<BaseModel> model(new Model_QuadrupoleLD("traj/quench"));
    int maxIter = 10000;
    int iter = 0;
    std::ofstream os;
    
    for (int i = 0; i < 100; i++) {

        model->createInitialState();
        iter = 0;
        while (iter < maxIter && !model->terminate()) {
            model->run(2);
            if ((iter + 1) % 100 == 0 || iter == 0) {
                state = model->getCurrState();
                os << iter << "\t" << state[0] << "\t" << state[1] << std::endl;
            }
            iter++;
        }

    }
}
void testQLearningMT(char* filename2){

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
    
    int num_threads = 2;
    std::vector<std::shared_ptr<BaseModel>> models;
    for (int i = 0; i < num_threads; i++){
        std::stringstream ss;
        ss << i;
        models.push_back(std::shared_ptr<BaseModel>(new Model_QuadrupoleLD("traj/thread" + ss.str() + "quench")));
    }
    
    
    RLSolver_2DTableMT rlSolver(models, 2, message3, n_rows, n_cols, dx1, dx2, minx1, minx2,num_threads);
//    rlSolver.getQTable().slice(2).ones();
    rlSolver.train();
}

void testQLearning(char* filename2){

    ReinforcementLearningParameter message2;
    QLearningSolverParameter message3;
    ReadProtoFromTextFile(filename2, &message2);
    message3 = message2.qlearningsolverparameter();
    std::shared_ptr<BaseModel> model(new Model_QuadrupoleLD("traj/"));
    
    int n_rows = 20;
    int n_cols = 20;
    
    double dx1 = 0.05;
    double dx2 = 0.05;
    double minx1 = 0;
    double minx2 = 0;
    
    RLSolver_2DTable rlSolver(model, 2, message3, n_rows, n_cols, dx1, dx2, minx1, minx2);
//    rlSolver.getQTable().slice(2).ones();
    rlSolver.train();
}