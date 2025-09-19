#include "Compartiments.h"

Arbres::Arbres() : Compartiments(), CI(370.23) {}

Arbres::Arbres(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille)
    : Compartiments(quantite, alpha, beta, gamma, delta, k, taille), CI(CI) {}

Arbres::Arbres(const Arbres& arbres) : Compartiments(arbres) {
    this->CI = arbres.CI;
}

Arbres::~Arbres() {}

double Arbres::S(double CT) const {
    // S(CT) = α * CT * (1 - CT/K)
    return getAlpha() * CT * (1.0 - CT / getK());
}

double Arbres::update(const Compartiments& arbre, const Compartiments& sol, const Compartiments& atmosphere, const Compartiments& humain, const Compartiments& ocean){
    // dCT/dt = S(CT) - β*CT - δ*CT - γ*CT

    double sequestration = S(arbre.getQuantite()) ;                 // séquestration par les arbres
    double respiration = getBeta() * arbre.getQuantite();         // respiration des arbres (β * CT)
    double transfertVersSol = getDelta() * arbre.getQuantite();   // transfert vers le sol (δ * CT)
    double litiere = getGamma() * arbre.getQuantite();            // litière (γ * CT)

    double dCT_dt = sequestration - respiration - transfertVersSol - litiere;
    return dCT_dt;
}