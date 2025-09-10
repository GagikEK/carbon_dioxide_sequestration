#include "atmosphere.h"
#include "arbres.h"

Atmosphere::Atmosphere() : Compartiments(), CI(0.0) {}

Atmosphere::Atmosphere(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille)
    : Compartiments(quantite, alpha, beta, gamma, delta, k, taille), CI(CI) {}

Atmosphere::~Atmosphere() {}

double Atmosphere::update(double CA, double CT, double CS) {
    // dCA/dt = -S(CT) + β*CT + δ*CS
    
    Arbres Arbre;
    double sequestration = Arbre.S(CT);       // carbone séquestré par les arbres
    double respirationArbres = getBeta() * CT;    // respiration des arbres (β * CT)
    double respirationSol = getDelta() * CS;      // respiration du sol (δ * CS)

    double dCA_dt = -sequestration + respirationArbres + respirationSol;
    return dCA_dt;
}