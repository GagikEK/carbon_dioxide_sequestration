#include "integration.h"

const double h = 1.0;  //temporaire le temps de créer la class "Simulation"

Eigen::VectorXd Integration::trapeze(Eigen::VectorXd fn, Eigen::VectorXd fn1){
    return h * (fn + fn1)/2; 
}

Eigen::VectorXd Integration::rectangle_droite(Eigen::VectorXd fn1){
    return h * fn1;
}

Eigen::VectorXd Integration::rectangle_gauche(Eigen::VectorXd fn){
    return h * fn;
}