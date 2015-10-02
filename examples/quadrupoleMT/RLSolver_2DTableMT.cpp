#include "RLSolver_2DTableMT.h"

using namespace ReinforcementLearning;


bool RLSolver_2DTableMT::finish_global = false;
int RLSolver_2DTableMT::threshFinishCount_global = 0;
std::mutex RLSolver_2DTableMT::QTable_mutex;


RLSolver_2DTableMT::RLSolver_2DTableMT(std::vector<std::shared_ptr<BaseModel>> m, int Dim, 
        ReinforcementLearning::QLearningSolverParameter para, int n_row0, int n_col0, 
        double dx, double dy, double min_x, double min_y, int num_threads0):
    RLSolver_2DTable(m[0],Dim,para,n_row0,n_col0,dx,dy,min_x,min_y),num_threads(num_threads0),models(m){
    
    

}
/*
 brief:
 * during the training, multiple threads will be initialized
 */
void RLSolver_2DTableMT::train() {
    int QTableOutputInterval = trainingPara.qtableoutputinterval();
    int ExperienceReplayInterval = trainingPara.experiencereplayinterval();
    experienceStopCriterion = trainingPara.experiencestopcriterion();
    std::thread *threads = new std::thread[num_threads];
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++){
        threads[thread_idx] = std::thread(RLSolver_2DTableMT::trainOnMT, models[thread_idx], thread_idx, this->trainingPara);    
    }
    
    while (RLSolver_2DTableMT::threshFinishCount_global < this->num_threads) {
        // after an episode, do experience reply
        if (RLSolver_2DTableMT::experienceVec.size() > 0 &&
                RLSolver_2DTableMT::experienceVec.size()%ExperienceReplayInterval == 0){
            std::cout << "Main Thread experience set size: "<< RLSolver_2DTableMT::experienceVec.size() << std::endl;
            this->replayExperience();
        }
        
        if (RLSolver_2DTableMT::experienceVec.size() > experienceStopCriterion) {
            finish_global = true;    
            break;
        }   
        
        if (RLSolver_2DTableMT::experienceVec.size() > 0 &&
                RLSolver_2DTableMT::experienceVec.size()%QTableOutputInterval == 0){
            for( int i = 0; i < RLSolver_2DTableMT::numActions; i++) {
                std::stringstream ss;        
                ss << (RLSolver_2DTableMT::experienceVec.size()/QTableOutputInterval);
                this->outputQ("QTable_" + ss.str() + "iter");
            }
        }
     }
 
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++){
        threads[thread_idx].join();    
    }
    
    std::cout << "experience set size: "<< RLSolver_2DTableMT::experienceVec.size() << std::endl;
        
        this->outputQ("QTableFinal");
        this->outputPolicy();
//        this->outputExperience("experience.dat");
        
}

void RLSolver_2DTableMT::trainOnMT(std::shared_ptr<BaseModel> m,int thread_idx,ReinforcementLearning::QLearningSolverParameter trainingPara){
    int iter;
    double maxQ, reward;
    int action;
    double epi = trainingPara.epsilon();
    int maxIter = trainingPara.numtrainingepisodes();
    int epiLength = trainingPara.episodelength();
    int controlInterval = trainingPara.controlinterval();
    int QTableOutputInterval = trainingPara.qtableoutputinterval();
    bool experienceReplayFlag = true;
    std::shared_ptr<BaseModel> model = m;
    std::shared_ptr<RandomStream> randChoice = std::make_shared<RandomStream>(0, model->getNumActions() - 1);
    for (int i = 0; i < maxIter; i++) {
        std::cout << "training Episodes " << i << " from thread " << thread_idx << std::endl;
        if (RLSolver_2DTableMT::finish_global){
            std::cout << "thread " << thread_idx << " stop!"<< std::endl;
            break;
        }
        iter = 0;
        model->createInitialState();
        while (!model->terminate() && iter < epiLength) {
            State oldState = model->getCurrState();
            if (randChoice->nextDou() < epi){
                RLSolver_2DTableMT::getMaxQ(oldState,&maxQ,&action);
            } else {
                action = randChoice->nextInt();
            }
            model->run(action, controlInterval);
            State newState = model->getCurrState();
            reward = model->getRewards();
            
            if (experienceReplayFlag) {
                // update the experienceVec with a lock
                std::unique_lock<std::mutex> lk(RLSolver_2DTableMT::QTable_mutex);
                RLSolver_2DTableMT::experienceVec.push_back(Experience(oldState, newState, action, reward));
                lk.unlock();
            }
             iter++;
        }
        
    }
    std::cout << "Thread " << thread_idx << " experience set size: "<< RLSolver_2DTableMT::experienceVec.size() << std::endl;
    std::unique_lock<std::mutex> lk(RLSolver_2DTableMT::QTable_mutex);
    RLSolver_2DTableMT::threshFinishCount_global +=1;
    lk.unlock();
           
}
