#ifndef Sol_H
#define Sol_H

#include <Eigen/Dense>
#include <iostream>
#include "compartiments.h"

class Sol : public Compartiments {
private:
    double CI;  // Initial condition

public:
    Sol();
    Sol(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille);
    ~Sol();

    double update(const Arbres& arbre, const Sol& sol, const Atmosphere& atmosphere, const Humains& humain, const Oceans& ocean) override;

    friend class Atmosphere;
    friend class Oceans;
    friend class Humains;
};

#endif