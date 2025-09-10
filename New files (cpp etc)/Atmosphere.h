#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include <Eigen/Dense>
#include <iostream>
#include "compartiments.h"

class Atmosphere : public Compartiments {
private:
    double CI;

public:
    Atmosphere();
    Atmosphere(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille);
    ~Atmosphere();

    double update(const Arbres& arbre, const Sol& sol, const Atmosphere& atmosphere, const Humains& humain, const Oceans& ocean) override;

    friend class Arbres;
    friend class Sol;
    friend class Oceans;
    friend class Humains;
};

#endif