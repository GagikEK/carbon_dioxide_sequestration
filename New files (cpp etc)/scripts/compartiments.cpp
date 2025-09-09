#include "compartiments.h"


double Compartiments::alpha  = 0.04;
double Compartiments::beta   = 0.02;
double Compartiments::_gamma = 0.015;
double Compartiments::delta  = 0.008;
double Compartiments::K      = 900;

Compartiments::Compartiments(){
    std::cout << "CREATION" << std::endl;
    this->quantite = 0;
    Compartiments::alpha = 0.04 ;
    Compartiments::beta = 0.02;
    Compartiments::_gamma = 0.015 ;
    Compartiments::delta = 0.008;
    Compartiments::K = 900;
    this->history = Eigen::VectorXd::Zero(100);

}

Compartiments::Compartiments(int quantite, double alpha, double beta, double gamma, double delta, double k, int taille){
    this->quantite = quantite;
    Compartiments::alpha = alpha;
    Compartiments::beta = beta;
    Compartiments::_gamma = gamma;
    Compartiments::delta = delta;
    Compartiments::K = k;

    this->history = Eigen::VectorXd::Zero(taille);

}

int Compartiments::getQuantite(){
    return this->quantite;
}

double Compartiments::getAlpha(){
    return Compartiments::alpha;
}

double Compartiments::getBeta(){
    return Compartiments::beta;
}

double Compartiments::getDelta(){
    return Compartiments::delta;
}

double Compartiments::getGamma(){
    return Compartiments::_gamma;
}


double Compartiments::getK(){
    return Compartiments::K;
}

Eigen::VectorXd Compartiments::getHistory(){
    return this->history;
}

void Compartiments::setConstantes(double alpha, double beta, double gamma, double delta, double K){
    Compartiments::alpha = alpha;
    Compartiments::beta = beta;
    Compartiments::_gamma = gamma;
    Compartiments::delta = delta;
    Compartiments::K = K;
}


void Compartiments::setQuantite(int quantite){
    this->quantite = quantite;
}

void Compartiments::setHistory(Eigen::VectorXd v){
    this->history = v;
}