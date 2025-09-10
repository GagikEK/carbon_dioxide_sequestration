#ifndef COMPARTIMENTS_H
#define COMPARTIMENTS_H

#include <Eigen/Dense>
#include <iostream>


class Compartiments{
    private :
        int quantite;
        static double alpha;
        static double beta;
        static double _gamma;
        static double delta;
        static double K;
        Eigen::VectorXd history;

    public:
        Compartiments();
        Compartiments(int quantite, double alpha, double beta, double gamma, double delta, double k, int taille);
        virtual ~Compartiments() = default;
        //virtual void update() = 0;
        //virtual void nextstep() = 0;

        int getQuantite() const;
        static double getAlpha();
        static double getBeta();
        static double getGamma();
        static double getDelta();
        static double getK();
        Eigen::VectorXd getHistory() const;

        static void setConstantes(double alpha, double beta, double gamma, double delta, double K);
        void setQuantite(int quantite);
        void setHistory(const Eigen::VectorXd& history);


        '''
        proposed by chat :
            // convenient helpers for history manipulation (no large replacements needed)
    void appendHistory(double value);      // push next value (increments history_pos)
    double lastValue() const;              // last pushed value (or 0.0 if empty)
    int historySize() const;               // capacity of history vector
    void reserveHistory(int size);         // allocate / resize history


        // (optional) an interface intended for polymorphism:
    // virtual void update() = 0;       // uncomment if you want derived classes to implement update()
    // virtual void nextstep() = 0;
        '''


    private:

        friend class Arbres;
        friend class Atmosphere;
        friend class Sol;
        friend class Oceans;
        friend class Humain;



};

#endif
