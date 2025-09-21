#ifndef COMPARTIMENTS_H
#define COMPARTIMENTS_H

#include <Eigen/Dense>
#include <iostream>

// forward declarations for compartment classes used in method signatures
class Arbres;
class Atmosphere;
class Sol;
class Oceans;
class Humains;

class Compartiments {
private:
    double quantite;
    static double alpha;
    static double beta;
    static double _gamma;
    static double delta;
    static double K;
    Eigen::VectorXd history;
    int history_pos;  // Current position in history

public:
    Compartiments();
    Compartiments(double quantite, double alpha, double beta, double gamma, double delta, double k, int taille);
    virtual ~Compartiments() = default;

    // Pure virtual method for updating compartment
    virtual double update(const Arbres& arbre, const Sol& sol, const Atmosphere& atmosphere, const Humains& humain, const Oceans& ocean) = 0;

    // Getters
    double getQuantite() const;
    static double getAlpha();
    static double getBeta();
    static double getGamma();
    static double getDelta();
    static double getK();
    Eigen::VectorXd getHistory() const;

    // Setters
    static void setConstantes(double alpha, double beta, double gamma, double delta, double K);
    void setQuantite(double quantite);
    void setHistory(const Eigen::VectorXd& history);

    // History management helpers
    void appendHistory(double value);
    double lastValue() const;
    int historySize() const;
    void reserveHistory(int size);
    void resetHistory();

private:
    friend class Arbres;
    friend class Atmosphere;
    friend class Sol;
    friend class Oceans;
    friend class Humains;
};

#endif