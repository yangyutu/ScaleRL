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
    int QTableOutputInterval = trainingPara.qtableoutputinterval(); //1000
    int ExperienceReplayInterval = trainingPara.experiencereplayinterval(); //30
    experienceStopCriterion = trainingPara.experiencestopcriterion();   //1e5最大experience次数
    int experienceReplayCounter = 0;
    int QTableOutputSizeCounter = 0;
    
    std::thread *threads = new std::thread[num_threads];
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++){
        threads[thread_idx] = std::thread(RLSolver_2DTableMT::trainOnMT, models[thread_idx], thread_idx, this->trainingPara);    
    }
    
    while (RLSolver_2DTableMT::threshFinishCount_global < this->num_threads) {  
        if (RLSolver_2DTableMT::experienceSetSize > experienceStopCriterion) { //Experience多于1e5次自动退出
            std::cout << "Main Thread: experience set satisfied: "<< RLSolver_2DTableMT::experienceSetSize << std::endl;            
            finish_global = true;    
            break;
        }   
        
        if (RLSolver_2DTableMT::experienceSetSize > QTableOutputSizeCounter){ //每存储1000次experience输出一次QTable
            QTableOutputSizeCounter += QTableOutputInterval;
            for( int i = 0; i < RLSolver_2DTableMT::numActions; i++) {
                std::stringstream ss;        
                ss << (RLSolver_2DTableMT::experienceSetSize/QTableOutputInterval);
                this->outputQ("QTable_" + ss.str() + "iter");
            }
            std::stringstream ss;        
            ss << (RLSolver_2DTableMT::experienceSetSize/QTableOutputInterval);
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
//  每个线程每次过程
    for (int i = 0; i < maxIter; i++) {
        std::cout << "training Episodes " << i << " from thread " << thread_idx << std::endl;
        if (RLSolver_2DTableMT::finish_global){
            std::cout << "thread " << thread_idx << " stop!"<< std::endl;
            break;
        }
        iter = 0;
//  1. 读取初始坐标状态
        model->createInitialState();
//  2. 判断结晶或到最大时间
        while (!model->terminate() && iter < epiLength) {
//      2.1 读取当前OP            
            State oldState = model->getCurrState(); 
//      2.2 决定OPT            
            if (randChoice->nextDou() < epi){ 
                std::unique_lock<std::mutex> lk(RLSolver_2DTableMT::QTable_mutex);
                RLSolver_2DTableMT::getMaxQ(oldState,&maxQ,&action); 
                lk.unlock();            
            } else { 
                action = randChoice->nextInt();
            }
//      2.3 模拟接下来的10s的运动，记录当前experience，更新QTable            
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
//      2.4 如果到时间需要进行ExperienceReplay，就重复记录并更新QTable                   
    std::cout << "Thread" << thread_idx << " replay experience set size: "<< 
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
