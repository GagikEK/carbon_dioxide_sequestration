#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <Eigen/Dense>

class Integration {
    public:
        static Eigen::VectorXd trapeze(Eigen::VectorXd fn, Eigen::VectorXd fn1);
        static Eigen::VectorXd rectangle_gauche(Eigen::VectorXd fn);
        static Eigen::VectorXd rectangle_droite(Eigen::VectorXd fn1);

};

#endif