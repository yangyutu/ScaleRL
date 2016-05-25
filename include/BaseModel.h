#pragma once
#include <vector>
namespace ReinforcementLearning {

    typedef std::vector<double> State;
/*
    std::ostream& operator<< (std::ostream& os, const State &S){
        for (auto &e:S){
            os << e;
        }    
    }
*/    
    class BaseModel {
    public:
        virtual ~BaseModel(){}
        virtual void run(int action) = 0;
	virtual void run(int action, int steps){
		for (int i = 0; i < steps; i++){
			run(action);
		}
	};
        virtual State getCurrState() {
            return currState;
        }
        virtual void createInitialState() = 0;
        virtual int getNumActions(){ return numActions;}
        virtual double getRewards() {}
        virtual bool terminate() {}
	virtual double callpsi6(){return 0;}
    protected:
        State currState, prevState;
        int numActions;
        int stateDim;
    };
    
    struct Experience{
        State oldState, newState;
        int action;
        double reward;
        Experience(State old0, State new0, int a0, double c0):
        oldState(old0),newState(new0), action(a0), reward(c0)
        {}
        friend std::ostream& operator<< (std::ostream &out, const Experience &exp){
            for (auto &e:exp.oldState){
                out << e << "\t";
            }    
            for (auto &e:exp.newState){
                out << e << "\t";
            }
            out<< exp.action << "\t";
            out << exp.reward << "\t";
        }
        
    };
    
}
