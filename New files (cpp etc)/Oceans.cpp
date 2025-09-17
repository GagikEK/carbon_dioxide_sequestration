#include "Oceans.h"
#include "Oceans.h"

Oceans::Oceans() : Compartiments(), CI(0.0) {}

Oceans::Oceans(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille)
    : Compartiments(quantite, alpha, beta, gamma, delta, k, taille), CI(CI) {}

Oceans::~Oceans() {}

double Oceans::S2(double Co) {
    // S(CT) = α * CT * (1 - CT/K)
    return alpha2 * Co * (1.0 - Co / k2);
}


double Oceans::update(const Arbres& arbre, const Sol& sol, const Atmosphere& atmosphere, const Humains& humain, const Oceans& ocean) override {
    // dCo/dt =  epsilon*CA - omega*CO

    double sequestration = S(ocean.getQuantite());    // le taux de séquestration du carbone dans l’océean
    double effet_resp = omega*ocean.getQuantite();          //effet de respiration de l’océan vers l’atmosphère
    double dCo_dt = sequestration -effet_resp;
    return dCo_dt;
}