#ifndef Simulation_H
#define Simulation_H

#include <Eigen/Dense>

class Simulation {
    private:
    
        int max_iter;
        double eps;
        double h;
        
        Eigen::VectorXd C0;
        double t0;
        double tf;

        Eigen::MatrixXd global_history;
        Eigen::VectorXd C;
        int current_pos;

        Eigen::VectorXd (*f)(double tn, Eigen::VectorXd Cn);

        int nb_pas();

    public:
        Simulation();
        Simulation(Eigen::VectorXd C0, double t0, double tf, double h, int max_iter=10000, double eps=1e-6);
        ~Simulation();


        int getMaxIter();
        void setMaxIter(int max_iter);

        double getTolerance();
        void setTolerance(double eps);

        int getPas();
        void setPas(double h);

        Eigen::VectorXd getCondInitial();
        void setCondInitial(Eigen::VectorXd C0);

        double getSimulBeginTime();
        void setSimulBeginTime(double t0);

        double getSimulTime();
        void setSimulTime(double t);

        double getSimulEndTime();
        void setSimulEndTime(double tf);

        Eigen::MatrixXd getSimulHistory();

        Eigen::VectorXd getSystemState();
        void setSystemState(Eigen::VectorXd C);

        void setSimulFunction(Eigen::VectorXd (*f)(double t, Eigen::VectorXd));

};

#endif