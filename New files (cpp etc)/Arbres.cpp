#include "arbres.h"

Arbres::Arbres() : Compartiments(), CI(0.0) {}

Arbres::Arbres(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille)
    : Compartiments(quantite, alpha, beta, gamma, delta, k, taille), CI(CI) {}

Arbres::~Arbres() {}

double Arbres::S(double CT) {
    // S(CT) = α * CT * (1 - CT/K)
    return getAlpha() * CT * (1.0 - CT / getK());
}

double Arbres::update(double CA, double CT, double CS) {
    // dCT/dt = S(CT) - β*CT - δ*CT - γ*CT
    
    double sequestration = S(CT);                 // séquestration par les arbres
    double respiration = getBeta() * CT;         // respiration des arbres (β * CT)
    double transfertVersSol = getDelta() * CT;   // transfert vers le sol (δ * CT)
    double litiere = getGamma() * CT;            // litière (γ * CT)

    double dCT_dt = sequestration - respiration - transfertVersSol - litiere;
    return dCT_dt;
}