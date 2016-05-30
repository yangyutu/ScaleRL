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
    int ExperienceReplayInterval = trainingPara.experiencereplayinterval();
    experienceStopCriterion = trainingPara.experiencestopcriterion();
    int experienceReplayCounter = 0;
    int QTableOutputSizeCounter = 0;
    
    std::thread *threads = new std::thread[num_threads];
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++){
        threads[thread_idx] = std::thread(RLSolver_2DTableMT::trainOnMT, models[thread_idx], thread_idx, this->trainingPara);    
    }
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++){
	threads[thread_idx].join();
    }
    this->outputQ("QTableFinal");
    this->outputPolicy();
    std::cout << "Training completed!" << std::endl;
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
    for (int i = 0; i < maxIter; i++) {
        std::cout << "Thread " << thread_idx+1 << " cycle " << i+1 << " start.." << std::endl;
        iter = 0;
        model->createInitialState();
        while (!model->terminate() && iter < epiLength) {           
            State oldState = model->getCurrState(); 
            if (randChoice->nextDou() < epi){ 
                std::unique_lock<std::mutex> lk(RLSolver_2DTableMT::QTable_mutex);
                RLSolver_2DTableMT::getMaxQ(oldState,&maxQ,&action); 
                lk.unlock();            
            } else { 
                action = randChoice->nextInt();
            }            
            model->run(action, controlInterval);
            State newState = model->getCurrState();
            reward = model->getRewards();
            std::unique_lock<std::mutex> lk(RLSolver_2DTableMT::QTable_mutex);
            RLSolver_2DTable::updateQ(Experience(oldState, newState, action, reward));
            RLSolver_2DTableMT::experienceVec.push_back(Experience(oldState, newState, action, reward));
            RLSolver_2DTableMT::experienceSetSize = RLSolver_2DTableMT::experienceVec.size();           
            lk.unlock();
            iter++;
        }                
    std::unique_lock<std::mutex> lk(RLSolver_2DTableMT::QTable_mutex);
    RLSolver_2DTableMT::replayExperience(RLSolver_2DTableMT::experienceSetSize);
    lk.unlock();
    std::cout << "Thread" << thread_idx+1 << " cycle "<< i+1 << " completed!"<< std::endl;
    }
}
