#include <atomic>
#include "RLSolver_2DTableMT.h"
using namespace ReinforcementLearning;


bool RLSolver_2DTableMT::finish_global = false;
std::atomic<int> RLSolver_2DTableMT::threshFinishCount_global(0);
std::mutex RLSolver_2DTableMT::QTable_mutex;
int RLSolver_2DTableMT::experienceSetSize = 0;

RLSolver_2DTableMT::RLSolver_2DTableMT(std::vector<std::shared_ptr<BaseModel>> m, int Dim, 
        ReinforcementLearning::QLearningSolverParameter para, int n_row0, int n_col0, 
        double dx, double dy, double min_x, double min_y, int num_threads0):
    RLSolver_2DTable(m[0],Dim,para,n_row0,n_col0,dx,dy,min_x,min_y),num_threads(num_threads0),models(m){}

void RLSolver_2DTableMT::train() {
    int QTableOutputInterval = trainingPara.qtableoutputinterval();
    // Record intermediate Q table every 1000 steps
    experienceStopCriterion = trainingPara.experiencestopcriterion();
    int QTableOutputSizeCounter = 0;
    std::thread *threads = new std::thread[num_threads];
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++){
        threads[thread_idx] = std::thread(RLSolver_2DTableMT::trainOnMT, models[thread_idx], thread_idx, this->trainingPara);    
    }
    
    while (RLSolver_2DTableMT::threshFinishCount_global < this->num_threads) {  
        if (RLSolver_2DTableMT::experienceSetSize > experienceStopCriterion) {
            std::cout << "Main Thread: experience set satisfied: "<< RLSolver_2DTableMT::experienceSetSize << std::endl;            
            finish_global = true;    
            break;
        }
     }
 
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++){threads[thread_idx].join();}
    
    std::cout << "Main thread experience set size: "<< RLSolver_2DTableMT::experienceVec.size() << std::endl;
        
    this->outputQ("QTableFinal");
    this->outputPolicy();
    std::cout << "Main thread: Training finish!" << std::endl;
    delete[] threads;   
}

void RLSolver_2DTableMT::trainOnMT(std::shared_ptr<BaseModel> m,int thread_idx,ReinforcementLearning::QLearningSolverParameter trainingPara){
    int iter;
    double maxQ, reward;
    int action;
    double epi = trainingPara.epsilon();
    int maxIter = trainingPara.numtrainingepisodes();
    int epiLength = trainingPara.episodelength();
    int controlInterval = trainingPara.controlinterval();
    int ExperienceReplayInterval = trainingPara.experiencereplayinterval();
    std::shared_ptr<BaseModel> model = m;
    std::shared_ptr<RandomStream> randChoice = std::make_shared<RandomStream>(0, model->getNumActions() - 1);
//  For each cycle of assembly (of a total of 30 cycles)
    for (int i = 0; i < maxIter; i++) {
        std::cout << "training Episodes " << i << " from thread " << thread_idx << std::endl;
        if (RLSolver_2DTableMT::finish_global){
            std::cout << "thread " << thread_idx << " stop!"<< std::endl;
            break;
        }
        iter = 0;
//  Read in initial conditions
        model->createInitialState();
//  Make progress if not crystallized or not over 1000s
        while (!model->terminate() && iter < epiLength) {
//      Read in current OP           
            State oldState = model->getCurrState(); 
//      Determine action based on current OP + 5% randomness          
            if (randChoice->nextDou() < epi){ 
                std::unique_lock<std::mutex> lk(RLSolver_2DTableMT::QTable_mutex);
                RLSolver_2DTableMT::getMaxQ(oldState,&maxQ,&action); 
                lk.unlock();            
            } else { 
                action = randChoice->nextInt();
            }
//      Simulate for next period, the outcome OP and reward of this period are recorded
            model->run(action, controlInterval);
            State newState = model->getCurrState();
            reward = model->getRewards();
//      Update the experiences list
            std::unique_lock<std::mutex> lk(RLSolver_2DTableMT::QTable_mutex);
            RLSolver_2DTable::updateQ(Experience(oldState, newState, action, reward));
            RLSolver_2DTableMT::experienceVec.push_back(Experience(oldState, newState, action, reward));
            RLSolver_2DTableMT::experienceSetSize = RLSolver_2DTableMT::experienceVec.size();           
            lk.unlock();
//      Move on to next period
            iter++;
        } 
// Replay the last 1000 experiences after each cycle
    std::cout << "Thread" << thread_idx << " has a total experiences of: "<< 
        RLSolver_2DTableMT::experienceSetSize << std::endl;
    std::unique_lock<std::mutex> lk(RLSolver_2DTableMT::QTable_mutex);
    RLSolver_2DTableMT::replayExperience(RLSolver_2DTableMT::experienceSetSize);
    lk.unlock();
    }
    std::cout << "Thread " << thread_idx << " FINISH with experience set size: "<< RLSolver_2DTableMT::experienceSetSize << std::endl;
    std::unique_lock<std::mutex> lk(RLSolver_2DTableMT::QTable_mutex);
    std::atomic_fetch_add(&RLSolver_2DTableMT::threshFinishCount_global,1);
    lk.unlock();
}

void RLSolver_2DTableMT::outputExperience(std::string filename) {
    
    std::ofstream os;
    os.open(filename);
    for (int i = 0; i < RLSolver_2DTableMT::experienceSetSize; i++){
        os << RLSolver_2DTableMT::experienceVec[i] << std::endl;
    }
 }
