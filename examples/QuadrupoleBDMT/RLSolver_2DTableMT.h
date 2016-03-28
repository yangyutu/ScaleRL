#pragma once
#include <utility>
#include <thread>
#include <mutex>
#include <atomic>
#include "common_RL.h"
#include "RLSolver_2DTable.h"


namespace ReinforcementLearning {

    class RLSolver_2DTableMT : public RLSolver_2DTable{
    public:
        RLSolver_2DTableMT(std::vector<std::shared_ptr<BaseModel>> m, int Dim, 
        ReinforcementLearning::QLearningSolverParameter para, int n_row0, int n_col0, 
        double dx, double dy, double min_x, double min_y, int num_threads);
        virtual ~RLSolver_2DTableMT() {}
        virtual void train();
    protected:
        void outputExperience(std::string filename);
        static void trainOnMT(std::shared_ptr<BaseModel> m, int idx, ReinforcementLearning::QLearningSolverParameter para);        
        static std::mutex QTable_mutex;
        std::vector<std::shared_ptr<BaseModel>> models;
        static bool finish_global;
        static std::atomic<int> threshFinishCount_global;
        int num_threads;        
        int experienceStopCriterion;
        static int experienceSetSize;
    };
}
