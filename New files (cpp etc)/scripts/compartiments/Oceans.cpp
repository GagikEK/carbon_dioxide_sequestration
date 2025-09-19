#include "Compartiments.h"

double omega = 0.0003;
double alpha2 = 0.001;
double k2 = 40000;


Oceans::Oceans() : Compartiments(), CI(40000) {}

Oceans::Oceans(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille)
    : Compartiments(quantite, alpha, beta, gamma, delta, k, taille), CI(CI) {}

Oceans::Oceans(const Oceans& oceans) : Compartiments(oceans){
    this->alpha2 = oceans.alpha2;
    this->k2 = oceans.k2;
    this->omega = oceans.omega;
    this->CI = oceans.CI;
}

Oceans::~Oceans() {}

double Oceans::S(double Co) {
    // S(CT) = α * CT * (1 - CT/K)
    return alpha2 * Co * (1.0 - Co / k2);
}


double Oceans::update(const Compartiments& arbre, const Compartiments& sol, const Compartiments& atmosphere, const Compartiments& humain, const Compartiments& ocean){
    // dCo/dt =  epsilon*CA - omega*CO

    double sequestration = S(ocean.getQuantite());    // le taux de séquestration du carbone dans l’océean
    double effet_resp = omega*ocean.getQuantite();          //effet de respiration de l’océan vers l’atmosphère
    double dCo_dt = sequestration -effet_resp;
    return dCo_dt;
}