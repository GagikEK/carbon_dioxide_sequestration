#ifndef ARBRES_H
#define ARBRES_H

#include <Eigen/Dense>
#include <iostream>
#include "compartiments.h"

class Arbres : public Compartiments {
private:
    double CI;  // Initial condition

public:
    Arbres();
    Arbres(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille);
    ~Arbres();

    double update(const Arbres& arbre, const Sol& sol, const Atmosphere& atmosphere, const Humains& humain, const Oceans& ocean) override;
    double S(double CT) const;  // Sequestration function

    friend class Atmosphere;
    friend class Sol;
    friend class Oceans;
    friend class Humains;
};

#endif