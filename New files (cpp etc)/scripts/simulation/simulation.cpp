#include "simulation.h"

Simulation::Simulation(){
    this->C0 = Eigen::VectorXd(4);
    this->C0 << 3306.296, //CA0 
                370.23,   //CT0  
                490.77,   //CS0
                40000;    //CO0

    this->t0 = 2024;
    this->tf = 2100;

    this->eps = 1e-6;
    this->max_iter = 10000;

    this->global_history = Eigen::MatrixXd::Zero(C0.size(), this->nb_pas()+1);
    this->global_history.col(0) = C0;
    this->current_pos = 1;

}

int Simulation::nb_pas(){
    int cpt = 0;
    int t = this->t0;
    int tf = this->tf;

    while(t < tf){
        cpt++;
        t = t+this->h;
    }
    return cpt;
}