#include "Compartiments.h"

Atmosphere::Atmosphere() : Compartiments(), CI(0.0) {}

Atmosphere::Atmosphere(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille)
    : Compartiments(quantite, alpha, beta, gamma, delta, k, taille), CI(CI) {}

Atmosphere::Atmosphere(const Atmosphere& atmosphere) : Compartiments(atmosphere) {
    this->CI = atmosphere.CI;
}

Atmosphere::~Atmosphere() {}

double Atmosphere::update(const Compartiments& arbre, const Compartiments& sol, const Compartiments& atmosphere, const Compartiments& humain, const Compartiments& ocean){
    // dCA/dt = -S(CT) + β*CT + δ*CS

    double sequestration = arbre.S(arbre.getQuantite());       // carbone séquestré par les arbres
    double respirationArbres = getBeta() * arbre.getQuantite();    // respiration des arbres (β * CT)
    double respirationSol = getDelta() * sol.getQuantite();      // respiration du sol (δ * CS)

    double dCA_dt = -sequestration + respirationArbres + respirationSol;
    return dCA_dt;
}