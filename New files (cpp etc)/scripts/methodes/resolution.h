#ifndef RESOLUTION_H
#define RESOLUTION_H

#include <Eigen/Dense>
#include <functional>

class Resolution{
    private:
        static double eps;
        static int max_iter;
    public:
        // use std::function callbacks so lambdas can be passed
        static Eigen::VectorXd pointFixe(const Eigen::VectorXd& C0, const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& F);
        static Eigen::VectorXd newton(const Eigen::VectorXd& C0, const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& F, const std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& dF);
        
        static double getEps();
        static int getMaxIter();

        static void setEps(double);
        static void setMaxIter(int);
};

#endif