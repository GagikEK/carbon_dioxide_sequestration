#define h 1 //temporaire le temps de créer la class "Simulation"

#include "integration.h"


Eigen::VectorXd Integration::trapeze(Eigen::VectorXd fn, Eigen::VectorXd fn1){
    return h * (fn + fn1)/2; 
}

Eigen::VectorXd Integration::rectangle_droite(Eigen::VectorXd fn1){
    return h * fn1;
}

Eigen::VectorXd Integration::rectangle_gauche(Eigen::VectorXd fn){
    return h * fn;
}