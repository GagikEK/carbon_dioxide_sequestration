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
        
        //virtual void update() = 0;
        //virtual void nextstep() = 0;

        int getQuantite();
        static double getAlpha();
        static double getBeta();
        static double getGamma();
        static double getDelta();
        static double getK();
        Eigen::VectorXd getHistory();

        static void setConstantes(double alpha, double beta, double gamma, double delta, double K);
        void setQuantite(int quantite);
        void setHistory(Eigen::VectorXd history);
};

#endif
