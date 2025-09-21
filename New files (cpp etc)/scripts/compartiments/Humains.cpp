#include "Humains.h"
#include "arbres.h"
#include <cmath>

// define static members
double Humains::a = 0.0;
double Humains::b = 0.0;

Humains::Humains() : Compartiments(), CI(0.0), t0(0.0), t(0.0) {}

Humains::Humains(double CI_, double t0_, double t_, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille)
    : Compartiments(quantite, alpha, beta, gamma, delta, k, taille), t0(t0_), t(t_), CI(CI_) {}

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

double Humains::update(const Arbres& arbre, const Sol& sol, const Atmosphere& atmosphere, const Humains& humain, const Oceans& ocean) {
    double CHT = std::exp(humain.getA()*humain.getT() + humain.getB())/1000000000.0;
    return CHT;
}