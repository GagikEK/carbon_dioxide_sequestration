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

    double update(double CA, double CT, double CS) override;
    double S(double CT);  // Sequestration function

    friend class Atmosphere;
    friend class Sol;
    friend class Oceans;
    friend class Humains;
};

#endif