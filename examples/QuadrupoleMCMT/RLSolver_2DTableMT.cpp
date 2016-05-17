#include <atomic>
#include "RLSolver_2DTableMT.h"
using namespace ReinforcementLearning;

std::atomic<int> RLSolver_2DTableMT::threshFinishCount_global(0);
std::mutex RLSolver_2DTableMT::QTable_mutex;
int RLSolver_2DTableMT::experienceSetSize = 0;

RLSolver_2DTableMT::RLSolver_2DTableMT(std::vector<std::shared_ptr<BaseModel>> m, int Dim, 
        ReinforcementLearning::QLearningSolverParameter para, int n_row0, int n_col0, 
        double dx, double dy, double min_x, double min_y, int num_threads0):
    RLSolver_2DTable(m[0],Dim,para,n_row0,n_col0,dx,dy,min_x,min_y),num_threads(num_threads0),models(m){}

void RLSolver_2DTableMT::train() {
    std::thread *threads = new std::thread[num_threads];
// Create threads and start simulations
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++){
        threads[thread_idx] = std::thread(RLSolver_2DTableMT::trainOnMT, models[thread_idx], thread_idx, this->trainingPara);    
    }
// After Simulation finished for all threads, join all results and yield QTable and policy
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++){
	threads[thread_idx].join();
    }
    this->outputQ("QTableFinal");
    this->outputPolicy();
    std::cout << "Training completed" << std::endl;
    delete[] threads;   
}

void RLSolver_2DTableMT::trainOnMT(std::shared_ptr<BaseModel> m,int thread_idx,ReinforcementLearning::QLearningSolverParameter trainingPara){
    int iter,action;
    double maxQ, reward;
    double epi = trainingPara.epsilon();
    int CycleNumber = trainingPara.numtrainingepisodes(); //number of cycles
    int CycleTimeLength = trainingPara.episodelength();
    int UpdateTime = trainingPara.controlinterval();
    std::shared_ptr<BaseModel> model = m;
    std::shared_ptr<RandomStream> randChoice = std::make_shared<RandomStream>(0, model->getNumActions() - 1);
/* For each cycle, run the simulation for a time length of CycleTimeLength (1000s); 
 * the control action is decided at the beginning (control point) of each control episode, 
 * which has a time length of UpdateTime (60s). The QTable is updated immediately after 
 * each episode. After each full cycle of simulation, the last 1000 experiences are 
 * recalled and used to update QTable again. 
*/
    for (int i = 0; i < CycleNumber; i++) {
        iter = 0;
        model->createInitialState();
// Termination criteria are either crystallization or 1000s of simulation
        while (!model->terminate() && iter < CycleTimeLength) {
// At the beginning of each control point, read in current state information and decide control action
            State oldState = model->getCurrState();            
            if (randChoice->nextDou() < epi){ 
                std::unique_lock<std::mutex> lk(RLSolver_2DTableMT::QTable_mutex);
                RLSolver_2DTableMT::getMaxQ(oldState,&maxQ,&action); 
                lk.unlock();            
            } else { 
                action = randChoice->nextInt();
            }
// Simulate through the control episode with a certain action, the end state is recorded to reflect the performance of the current action using at the start state; the QTable is updated immediately for this experience.
            model->run(action, UpdateTime);
            State newState = model->getCurrState();
            reward = model->getRewards();
            std::unique_lock<std::mutex> lk(RLSolver_2DTableMT::QTable_mutex);
            RLSolver_2DTable::updateQ(Experience(oldState, newState, action, reward));
            RLSolver_2DTableMT::experienceVec.push_back(Experience(oldState, newState, action, reward)); 
            lk.unlock();
            iter++;
        } 
// After the full cycle, the last 1000 experiences are recalled and used to update QTable again
    std::unique_lock<std::mutex> lk(RLSolver_2DTableMT::QTable_mutex);
    RLSolver_2DTableMT::experienceSetSize = RLSolver_2DTableMT::experienceVec.size();
    RLSolver_2DTableMT::replayExperience(RLSolver_2DTableMT::experienceSetSize);
    lk.unlock();
    std::cout << "Thread " << thread_idx+1 << " cycle " << i+1 << " completed" << std::endl;
    }
    std::cout << "Thread " << thread_idx+1 << " completed" << std::endl;
    std::unique_lock<std::mutex> lk(RLSolver_2DTableMT::QTable_mutex);
    std::atomic_fetch_add(&RLSolver_2DTableMT::threshFinishCount_global,1);
    lk.unlock();
           
}
