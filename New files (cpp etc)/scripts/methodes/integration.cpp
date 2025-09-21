#include "integration.h"

Eigen::VectorXd Integration::trapeze(const Eigen::VectorXd& fn, const Eigen::VectorXd& fn1, double h){
    return h * (fn + fn1) / 2.0;
}

Eigen::VectorXd Integration::rectangle_droite(const Eigen::VectorXd& fn1, double h){
    return h * fn1;
}

Eigen::VectorXd Integration::rectangle_gauche(const Eigen::VectorXd& fn, double h){
    return h * fn;
}