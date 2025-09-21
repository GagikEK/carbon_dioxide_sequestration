#ifndef COMPARTIMENTS_H
#define COMPARTIMENTS_H

#include <Eigen/Dense>
#include <cmath>

class Compartiments {
protected:
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
    Compartiments(const Compartiments& compartiments);
    virtual ~Compartiments() = default;

    // Pure virtual method for updating compartment
    virtual double update(const Compartiments& arbre, const Compartiments& sol, const Compartiments& atmosphere, const Compartiments& humain, const Compartiments& ocean) = 0;
    virtual double S(double quantite) const = 0 ; 

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



class Atmosphere : public Compartiments {
private:
    double CI;

public:
    Atmosphere();
    Atmosphere(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille);
    Atmosphere(const Atmosphere& atmosphere);

    ~Atmosphere();

    double update(const Compartiments& arbre, const Compartiments& sol, const Compartiments& atmosphere, const Compartiments& humain, const Compartiments& ocean) override;
    double S(double CT) const override;

    friend class Arbres;
    friend class Sol;
    friend class Oceans;
    friend class Humains;
};



class Arbres : public Compartiments {
private:
    double CI;  // Initial condition

public:
    Arbres();
    Arbres(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille);
    Arbres(const Arbres& arbre);
    
    ~Arbres();

    double update(const Compartiments& arbre, const Compartiments& sol, const Compartiments& atmosphere, const Compartiments& humain, const Compartiments& ocean) override;
    double S(double CT) const override;

    friend class Atmosphere;
    friend class Sol;
    friend class Oceans;
    friend class Humains;
};



class Sol : public Compartiments {
private:
    double CI;  // Initial condition

public:
    Sol();
    Sol(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille);
    Sol(const Sol& sol);
    
    ~Sol();

    double update(const Compartiments& arbre, const Compartiments& sol, const Compartiments& atmosphere, const Compartiments& humain, const Compartiments& ocean) override;
    double S(double CT) const override;

    friend class Atmosphere;
    friend class Sol;
    friend class Oceans;
    friend class Humains;
};



class Oceans : public Compartiments {
protected:
    static double omega;
    static double alpha2;
    static double k2;
    double CI;

public:
    Oceans();
    Oceans(double CI, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille);
    Oceans(const Oceans& oceans);

    ~Oceans();

    double update(const Compartiments& arbre, const Compartiments& sol, const Compartiments& atmosphere, const Compartiments& humain, const Compartiments& ocean) override;
    double S(double CT) const override;

    friend class Arbres;
    friend class Sol;
    friend class Oceans;
    friend class Humains;
};



class Humains : public Compartiments {
private:
    static double a;
    static double b;
    double t0;
    double t;
    double CI;

public:
    Humains();
    Humains(double CI,double t0, double t, double quantite, double alpha, double beta, double gamma, double delta, double k, int taille);
    Humains(const Humains& humains);
    
    ~Humains();

    double update(const Compartiments& arbre, const Compartiments& sol, const Compartiments& atmosphere, const Compartiments& humain, const Compartiments& ocean) override;
    double S(double CT) const override;

    double getT0() const;
    double getT() const;
    static double getA();
    static double getB();

    friend class Arbres;
    friend class Sol;
    friend class Oceans;
    friend class Humains;
};
#endif