#ifndef RESOLUTION_H
#define RESOLUTION_H

#include <Eigen/Dense>

class Resolution{
    private:
        static double eps;
        static int max_iter;
    public:
        static Eigen::VectorXd pointFixe(Eigen::VectorXd C0, Eigen::VectorXd (*F)(Eigen::VectorXd));
        static Eigen::VectorXd newton(Eigen::VectorXd C0, Eigen::VectorXd (*F)(Eigen::VectorXd), Eigen::MatrixXd (*dF)(Eigen::VectorXd));
        
        static double getEps();
        static int getMaxIter();

        static void setEps(double);
        static void setMaxIter(int);
};

#endif