#include "Compartiments.h"

// Initialize static members
double Compartiments::alpha = 0.04;
double Compartiments::beta = 0.02;
double Compartiments::_gamma = 0.015;
double Compartiments::delta = 0.008;
double Compartiments::K = 900;

Compartiments::Compartiments() : quantite(0.0), history_pos(0) {
    this->history = Eigen::VectorXd::Zero(100);
}

Compartiments::Compartiments(double quantite, double alpha, double beta, double gamma, double delta, double k, int taille)
    : quantite(quantite), history_pos(0) {
    Compartiments::alpha = alpha;
    Compartiments::beta = beta;
    Compartiments::_gamma = gamma;
    Compartiments::delta = delta;
    Compartiments::K = k;
    this->history = Eigen::VectorXd::Zero(taille);
}

Compartiments::Compartiments(const Compartiments& compartiments){
    this->alpha = compartiments.getAlpha();
    this->beta = compartiments.getBeta();
    this->_gamma = compartiments.getGamma();
    this->delta = compartiments.getDelta();
    this->K = compartiments.getK();
    this->history = compartiments.getHistory();
    this->history_pos = compartiments.history_pos;
    this->quantite = compartiments.getQuantite();
}


double Compartiments::getQuantite() const {
    return this->quantite;
}

double Compartiments::getAlpha() {
    return Compartiments::alpha;
}

double Compartiments::getBeta() {
    return Compartiments::beta;
}

double Compartiments::getDelta() {
    return Compartiments::delta;
}

double Compartiments::getGamma() {
    return Compartiments::_gamma;
}

double Compartiments::getK() {
    return Compartiments::K;
}

Eigen::VectorXd Compartiments::getHistory() const {
    return this->history;
}

void Compartiments::setConstantes(double alpha, double beta, double gamma, double delta, double K) {
    Compartiments::alpha = alpha;
    Compartiments::beta = beta;
    Compartiments::_gamma = gamma;
    Compartiments::delta = delta;
    Compartiments::K = K;
}

void Compartiments::setQuantite(double quantite) {
    this->quantite = quantite;
}

void Compartiments::setHistory(const Eigen::VectorXd& v) {
    this->history = v;
}

void Compartiments::appendHistory(double value) {
    if (history_pos < history.size()) {
        history(history_pos) = value;
        history_pos++;
    }
}

double Compartiments::lastValue() const {
    if (history_pos > 0) {
        return history(history_pos - 1);
    }
    return 0.0;
}

int Compartiments::historySize() const {
    return history.size();
}

void Compartiments::reserveHistory(int size) {
    history = Eigen::VectorXd::Zero(size);
    history_pos = 0;
}

void Compartiments::resetHistory() {
    history_pos = 0;
    history.setZero();
}