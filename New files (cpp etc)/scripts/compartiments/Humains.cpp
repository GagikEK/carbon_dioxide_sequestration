#include "Compartiments.h"

double Humains::a = 0.022863115449118603;
double Humains::b = -21.769830005922582;

Humains::Humains() : Compartiments() {
    this->t0 = 2024; 
    this->t = 2024; 
    this->CI = std::exp(this->getA()*this->getT0() + this->getB())/1000000000;
    this->quantite = this->CI;
}

Humains::Humains(double CI,double t0, double t, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille)
    : Compartiments(quantite, alpha, beta, gamma, delta, k, taille), t0(t0),t(t), CI(CI) {}

Humains::Humains(const Humains& humains) : Compartiments(humains) {
    this->t0 = humains.t0;
    this->t = humains.t;
    this->a = humains.getA();
    this->b = humains.getB();
    this->CI = humains.CI;
}

Humains::~Humains() {}

double Humains::getA(){
    return a;
}
double Humains::getB(){
    return b;
}
double Humains::getT() const{
    return t;
}
double Humains::getT0() const{
    return t0;
}
double Humains::update(const Compartiments& arbre, const Compartiments& sol, const Compartiments& atmosphere, const Compartiments& humain, const Compartiments& ocean){
 
    double CHT = exp(this->getA()*this->getT() + this->getB())/1000000000;
    return CHT;
}

double Humains::S(double CT) const{
    return 0;
}