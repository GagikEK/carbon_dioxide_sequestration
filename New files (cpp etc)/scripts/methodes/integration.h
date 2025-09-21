#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <Eigen/Dense>

class Integration {
    public:
    // integration rules take function values and timestep h
    static Eigen::VectorXd trapeze(const Eigen::VectorXd& fn, const Eigen::VectorXd& fn1, double h);
    static Eigen::VectorXd rectangle_gauche(const Eigen::VectorXd& fn, double h);
    static Eigen::VectorXd rectangle_droite(const Eigen::VectorXd& fn1, double h);

};

#endif