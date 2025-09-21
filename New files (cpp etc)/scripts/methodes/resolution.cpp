#include "resolution.h"
#include <iostream>

int Resolution::max_iter = 10000;
double Resolution::eps = 1e-6;

double Resolution::getEps(){
    return Resolution::eps;
}

int Resolution::getMaxIter(){
    return Resolution::max_iter;
}

void Resolution::setEps(double eps){
    Resolution::eps = eps;
}

void Resolution::setMaxIter(int maxiter){
    Resolution::max_iter = maxiter;
}

Eigen::VectorXd Resolution::pointFixe(const Eigen::VectorXd& C0, const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& F){
    Eigen::VectorXd Cnk = C0;
    Eigen::VectorXd Cnk1 = C0;
    for(int i=0 ; i< Resolution::max_iter ; i++){
        Cnk1 = F(Cnk1);
        if((Cnk1 - Cnk).norm() < Resolution::eps){
            return Cnk1;
        }
        Cnk = Cnk1;
    }

    std::cerr << "Point fixe : non convergence après " << Resolution::max_iter << " itérations" << std::endl;

    return Cnk1;
}

Eigen::VectorXd Resolution::newton(const Eigen::VectorXd& C0, const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& F, const std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& dF){
    Eigen::VectorXd Cnk = C0;
    for(int i=0 ; i < Resolution::max_iter ; i++){
        if(F(Cnk).norm() < Resolution::eps){
            return Cnk;
        }
        Cnk = Cnk - dF(Cnk).inverse() * F(Cnk);
    }

    std::cerr << "Newton : non convergence après " << Resolution::max_iter << " itérations" << std::endl;

    return Cnk;
}