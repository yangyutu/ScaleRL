#pragma once
#include <utility>
#include "common_RL.h"
#include "RLSolverBase.h"


namespace ReinforcementLearning {

    class RLSolver_2DTable : public RLSolverBase{
    public:
        RLSolver_2DTable(std::shared_ptr<BaseModel> m, int Dim, 
        ReinforcementLearning::QLearningSolverParameter para, int n_row0, int n_col0, 
        double dx, double dy, double min_x, double min_y);
        virtual ~RLSolver_2DTable() {}
        virtual void train();
        virtual void test();
	static arma::cube CountMap;
        static void replayExperience();
        static void updateQ(Experience);
        static void getMaxQ(const State& S, double* Q, int* action);
        virtual arma::cube& getQTable(){return QTable;}
        virtual void loadQTable(std::string filetag);
    protected:
        void outputPolicy() const;
	void outputQ(std::string filename) const;
        void writeTrajectory(int iter, std::ostream &os, int action, State state, double reward) const;
        static std::pair<int, int> stateToIndex(const State & S);
        static arma::cube QTable;
        static int n_rows, n_cols, numActions;
        static double dx1, dx2, minx1, minx2;
        static arma::Mat<int> count;
	
        static std::vector<Experience> experienceVec;

    };
}
