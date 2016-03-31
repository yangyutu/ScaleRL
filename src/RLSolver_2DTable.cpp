#include "RLSolver_2DTable.h"

using namespace ReinforcementLearning;

arma::cube RLSolver_2DTable::QTable;
int RLSolver_2DTable::n_rows, RLSolver_2DTable::n_cols, RLSolver_2DTable::numActions;
double RLSolver_2DTable::dx1, RLSolver_2DTable::dx2, RLSolver_2DTable::minx1, RLSolver_2DTable::minx2;
arma::Mat<int> RLSolver_2DTable::count;
std::vector<Experience> RLSolver_2DTable::experienceVec;
arma::cube RLSolver_2DTable::CountMap;


RLSolver_2DTable::RLSolver_2DTable(std::shared_ptr<BaseModel> m, int Dim, 
        ReinforcementLearning::QLearningSolverParameter para, int n_row0, int n_col0, 
        double dx, double dy, double min_x, double min_y):
    RLSolverBase(m,Dim,para){
    
    n_rows = n_row0;
    n_cols = n_col0;    
    dx1 = dx;
    dx2 = dy;
    minx1 = min_x;
    minx2 = min_y;
    numActions = model->getNumActions();
    QTable.zeros(n_rows, n_cols, numActions);
    count.zeros(n_rows, n_cols);
    CountMap.zeros(n_rows, n_cols, 4);
}


void RLSolver_2DTable::train() {
    int iter;
    double maxQ, reward;
    int action;
    double epi = trainingPara.epsilon();
    int maxIter = trainingPara.numtrainingepisodes();
    int epiLength = trainingPara.episodelength();
    int controlInterval = trainingPara.controlinterval();
    int QTableOutputInterval = trainingPara.qtableoutputinterval();
    bool experienceReplayFlag = true;
    for (int i = 0; i < maxIter; i++) {
        std::cout << "training Episodes " << i << std::endl;
        iter = 0;
        model->createInitialState();
        while (!model->terminate() && iter < epiLength) {
            State oldState = model->getCurrState();
            if (randChoice->nextDou() < epi){
                getMaxQ(oldState,&maxQ,&action);
            } else {
                action = randChoice->nextInt();
            }
            model->run(action, controlInterval);
            State newState = model->getCurrState();
            reward = model->getRewards();
            this->updateQ(Experience(oldState, newState, action, reward));
            if (experienceReplayFlag) {
                experienceVec.push_back(Experience(oldState, newState, action, reward));
            }
             iter++;
        }
        // after an episode, do experience reply
        this->replayExperience();
        std::cout << "duration: " << iter << std::endl;
        if (i == 0 || ((i + 1) % QTableOutputInterval == 0)) {
            std::stringstream ss;        
            ss << i;
            this->outputQ("QTable_" + ss.str() + "iter");
        }
    }
        this->outputQ("QTableFinal");
        this->outputPolicy();
}

void RLSolver_2DTable::replayExperience(){
    // replay in reverse order
    for (auto exp = RLSolver_2DTable::experienceVec.rbegin(); exp != RLSolver_2DTable::experienceVec.rend(); ++exp){
        updateQ(*exp);
    }
}

void RLSolver_2DTable::test(){
// do some test
    int maxIter = trainingPara.numtrainingepisodes();
    int epiLength = trainingPara.episodelength();
    int controlInterval = trainingPara.controlinterval();
    int iter;
    int action;
    double maxQ;
    for (int i = 0; i < 0.1*maxIter; i++) {
        std::cout << "testing Episodes " << i << std::endl;
        iter = 0;
        model->createInitialState();
        while (!model->terminate() && iter < epiLength) {
            State oldState = model->getCurrState();          
            getMaxQ(oldState,&maxQ,&action);
            model->run(action, controlInterval);
            State newState = model->getCurrState();
            double reward = model->getRewards();
            iter++;
        }  
         std::cout << "duration: " << iter << std::endl;
    }   
}

void RLSolver_2DTable::updateQ(Experience exp) {
    double learningRate = trainingPara.learningrate();
    int action;
    double maxQ;
    double discount = trainingPara.discount();
    getMaxQ(exp.newState, &maxQ, &action);
    std::pair<int, int> index0 = stateToIndex(exp.oldState);
    count(index0.first,index0.second) += 1;
    CountMap(index0.first,index0.second,exp.action+1) += 1;
    QTable(index0.first,index0.second,exp.action) += 
            (1 / (1 + CountMap(index0.first,index0.second,exp.action))) * (exp.reward + discount * maxQ - QTable(index0.first, index0.second, exp.action));
}

void RLSolver_2DTable::getMaxQ(const State& S, double* maxQ, int* action){
    std::pair<int, int> index = RLSolver_2DTable::stateToIndex(S);
    double max = -std::numeric_limits<double>::max();
    for (int i = 0; i < numActions; i++){
        if(max < QTable(index.first,index.second, i)){
            max = QTable(index.first,index.second, i);
            *action = i;
        }
    }
    *maxQ = max;
}

std::pair<int, int> RLSolver_2DTable::stateToIndex(const State& S){
    int idx1, idx2;
    idx1 = (int) ((S[0] - minx1)/dx1);
    idx2 = (int) ((S[1] - minx2)/dx2);
    if (idx1 < 0) idx1 = 0;
    if (idx1 >= n_rows) idx1 = n_rows - 1; 
    if (idx2 < 0) idx2 = 0;
    if (idx2 >= n_cols) idx2 = n_cols - 1; 
    
    
    return std::pair<int, int>(idx1,idx2);
}

void RLSolver_2DTable::outputQ(std::string filename) const {
    
    for (int i = 0; i < numActions; ++i) {
        arma::mat temp = QTable.slice(i);
        std::stringstream ss;
        ss << i;
        temp.save(filename + ss.str() + ".dat", arma::raw_ascii);
    }
}



void RLSolver_2DTable::loadQTable(std::string filetag){
    for (int i = 0; i < numActions; ++i) {
        arma::mat temp;
        std::stringstream ss;
        ss << i;
        temp.load(filetag + ss.str() + ".dat");
	QTable.slice(i) = temp;    
	}
}

void RLSolver_2DTable::outputPolicy() const{
    arma::Mat<int> actionMap(n_rows, n_cols, arma::fill::zeros);
    arma::mat QMap(n_rows, n_cols, arma::fill::ones);
    QMap *= -1;
    double maxQ;
    int action;
    for (int i = 0; i < n_rows; i++){
        for (int j = 0; j < n_cols; j++){
            maxQ = -1000;
            action = 0;
            for (int m = 0; m < numActions; m++){
                if(maxQ < QTable(i, j, m)){
                    maxQ = QTable(i, j, m);
                    action = m;
                }
            }
            actionMap(i,j) = action;
            QMap(i,j) = maxQ;
            if (count(i,j) == 0){
                actionMap(i,j) = -1;
            }
        }
    }
    actionMap.save("actionMap.dat", arma::raw_ascii);
    QMap.save("QMap.dat", arma::raw_ascii);
    count.save("count.dat", arma::raw_ascii);
}

