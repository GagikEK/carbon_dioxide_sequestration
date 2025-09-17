#ifndef Humains_H
#define Humains_H

#include <Eigen/Dense>
#include <iostream>
#include "compartiments.h"

class Humains : public Compartiments {
private:
   static double a;
    static double b;
    double t0;
    double t;
    double CI;

public:
    Humains();
    Humains(double CI,double t0, double t, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille);
    ~Humains();

    double update(const Arbres& arbre, const Sol& sol, const Atmosphere& atmosphere, const Humains& humain, const Oceans& ocean) override;

    double getT0();
    double getT();
    static double getA();
    static double getB();

    friend class Arbres;
    friend class Sol;
    friend class Oceans;
    friend class Humains;
};

#endif