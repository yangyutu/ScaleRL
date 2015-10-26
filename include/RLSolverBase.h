#pragma once
#include "common_RL.h"
#include "BaseModel.h"
namespace ReinforcementLearning {

    class RLSolverBase {
    public:

        RLSolverBase(std::shared_ptr<BaseModel> m, int Dim, ReinforcementLearning::QLearningSolverParameter para) {
            trainingPara = para;
            model = m;
            stateDim = Dim;
            randChoice = std::make_shared<RandomStream>(0, model->getNumActions() - 1);
        }

        virtual ~RLSolverBase() {
        }
        virtual void train() = 0;
        virtual double getRewards(const State& newS) const {return 0.0;};
    protected:
        int stateDim;
        std::shared_ptr<BaseModel> model;
        std::shared_ptr<RandomStream> randChoice;
        static QLearningSolverParameter trainingPara;

    };
}
