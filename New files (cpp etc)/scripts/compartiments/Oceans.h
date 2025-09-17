#ifndef Oceans_H
#define Oceans_H

#include <Eigen/Dense>
#include <iostream>
#include "compartiments.h"

class Oceans : public Compartiments {
protected:
    static double omega=0.0003;
    static double alpha2=0.001;
    static double k2=40000;
    double CI;

public:
    Oceans();
    Oceans(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille);
    ~Oceans();

    double update(const Arbres& arbre, const Sol& sol, const Atmosphere& atmosphere, const Humains& humain, const Oceans& ocean) override;
    double S2(double CT); 

    friend class Arbres;
    friend class Sol;
    friend class Oceans;
    friend class Humains;
};

#endif