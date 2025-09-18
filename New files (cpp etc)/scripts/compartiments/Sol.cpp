#include "Compartiments.h"

Sol::Sol() : Compartiments(), CI(0.0) {}

Sol::Sol(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille)
    : Compartiments(quantite, alpha, beta, gamma, delta, k, taille), CI(CI) {}

Sol::Sol(const Sol& sol) : Compartiments(sol){
    this->CI = sol.CI;
}

Sol::~Sol() {}

double Sol::update(const Arbres& arbre, const Sol& sol, const Atmosphere& atmosphere, const Humains& humain, const Oceans& ocean) override {
    // dCs/dt =  _gamma*CT - delta*CS + delta*CT

    double litière = getGamma()*arbre.getQuantite();     // litière (γ * CT)
    double respirationSol = getDelta() * sol.getQuantite();      // respiration du sol (δ * CS)  
    double transfertVersSol = getDelta() * arbre.getQuantite();    // transfertVersSol (δ * CT)
    
    double dCS_dt = litière + respirationSol + transfertVersSol;
    return dCS_dt;
}