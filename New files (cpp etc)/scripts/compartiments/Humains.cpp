#include "Humains.h"
#include "arbres.h"
#include <cmath> 

Humains::Humains() : Compartiments(), CI(0.0),to(0.0),t(0.0) {}

Humains::Humains(double CI,double t0, double t, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille)
    : Compartiments(quantite, alpha, beta, gamma, delta, k, taille), to(to),t(t) {}

Humains::~Humains() {}

static double Humains::getA(){
    return a;
}
static double Humains::getB(){
    return B;
}
double Humains::getT(){
    return t;
}
double Humains::getT0(){
    return t0;
}
double Humains::update(const Arbres& arbre, const Sol& sol, const Atmosphere& atmosphere, const Humains& humain, const Oceans& ocean) override {
 
    double CHT = exp(humain.getA()*humain.getT() + humain.getB())/1000000000;
    return CHT;
}