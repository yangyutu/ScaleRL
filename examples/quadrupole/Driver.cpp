#include <iostream>
#include "boost/multi_array.hpp"
#include "common_RL.h"
#include "Model_QuadrupoleLD.h"
#include "RLSolver_2DTable.h"
using namespace ReinforcementLearning;

void testModel();
void testQLearning(char* filename);
int main(int argc, char* argv[]){
    testModel();
    testQLearning(argv[1]);
    return 1;
}


void testModel() {
    State state;
    std::shared_ptr<BaseModel> model(new Model_QuadrupoleLD);
    int maxIter = 10000;
    int iter = 0;
    std::ofstream os;
    
    for (int i = 0; i < 100; i++) {
        std::stringstream ss;
        ss << i;
        os.open("traj/quench"+ss.str()+".dat");
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
        os.close();
    }
}
void testQLearning(char* filename2){

    QLearningSolverParameter message3;
    ReadProtoFromTextFile(filename2, &message3);
    std::shared_ptr<BaseModel> model(new Model_QuadrupoleLD());
    
    int n_rows = 20;
    int n_cols = 78;
    
    double dx1 = 0.05;
    double dx2 = -100.0;
    double minx1 = 0;
    double minx2 = 24000;
    
    RLSolver_2DTable rlSolver(model, 2, message3, n_rows, n_cols, dx1, dx2, minx1, minx2);
//    rlSolver.getQTable().slice(2).ones();
    rlSolver.train();
}

