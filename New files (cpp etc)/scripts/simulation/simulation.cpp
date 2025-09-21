#include "simulation.h"
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <iomanip>
// methods
#include "../methodes/integration.h"
#include "../methodes/resolution.h"

// We'll provide a default RHS that maps the current compartment quantities to the
// ODE right-hand side vector [dCA/dt, dCT/dt, dCS/dt, dCO/dt].
// The Simulation object keeps compartment instances; here we provide a helper
// that builds the RHS given explicit compartment quantities. The Simulation
// class will construct a std::function wrapper that captures its compartments
// and calls this helper.
Eigen::VectorXd f_base(double tn, const Eigen::VectorXd& Cn) {
    // If called standalone, return zeros for safety
    return Eigen::VectorXd::Zero(Cn.size());
}

Simulation::Simulation()
    : arbre(), ocean(), humain(), sol(), atmosphere(),
      max_iter(10000), eps(1e-6), h(1.0), C0(), t0(0.0), tf(0.0),
      global_history(), C(), current_pos(0), time_vector(), f(f_base), integrator_(Integrator::rectangle_gauche)
{
    // Default f captures the compartments by reference and builds the 4-component RHS
    f = [this](double tn, const Eigen::VectorXd& Cn)->Eigen::VectorXd {
        // Build the RHS by delegating to each compartment's update() method.
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(4);

        // Each update(...) expects references to the other compartments so pass
        // the current objects. Their internal getQuantite() should reflect the
        // current state; ensure Simulation sets them before calling run.
        rhs(0) = atmosphere.update(arbre, sol, atmosphere, humain, ocean);
        rhs(1) = arbre.update(arbre, sol, atmosphere, humain, ocean);
        rhs(2) = sol.update(arbre, sol, atmosphere, humain, ocean);
        rhs(3) = ocean.update(arbre, sol, atmosphere, humain, ocean);

        return rhs;
    };
}

Simulation::Simulation(const Eigen::VectorXd& C0_, double t0_, double tf_, double h_, int max_iter_, double eps_, const Arbres& arbre_, const Sol& sol_, const Atmosphere& atmosphere_, const Humains& humain_, const Oceans& ocean_)
    : arbre(arbre_), ocean(ocean_), humain(humain_), sol(sol_), atmosphere(atmosphere_),
      max_iter(max_iter_), eps(eps_), h(h_), C0(C0_), t0(t0_), tf(tf_),
    global_history(), C(C0_), current_pos(0), time_vector(), f(f_base), integrator_(Integrator::rectangle_gauche)
{
    if (h <= 0) throw std::invalid_argument("step size h must be > 0");
    if (tf < t0) throw std::invalid_argument("end time must be >= start time");

    int rows = nb_pas() + 1;
    global_history = Eigen::MatrixXd::Zero(std::max(1, rows), std::max(1, (int)C.size()));
    if (C.size() > 0) global_history.row(0) = C.transpose();

    time_vector = Eigen::VectorXd::Zero(rows + 0);
    if (rows > 0) {
        for (int i = 0; i <= nb_pas(); ++i) time_vector(i) = t0 + i * h;
    }
}

Simulation::~Simulation() = default;

int Simulation::nb_pas() const {
    if (h <= 0) return 0;
    return static_cast<int>(std::floor((tf - t0) / h + 0.5));
}

void Simulation::setSimulationMethod(const std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)>& func) {
    f = func;
}

void Simulation::setSimulFunction(Eigen::VectorXd (*fp)(double, const Eigen::VectorXd&)) {
    if (fp) f = fp; else f = f_base;
}

void Simulation::setIntegrator(Integrator integrator_) { integrator_ = integrator_; }
Simulation::Integrator Simulation::getIntegrator() const { return integrator_; }

void Simulation::run() {
    int steps = nb_pas();
    if (steps <= 0) return;

    double t = t0;
    global_history = Eigen::MatrixXd::Zero(steps + 1, C.size());
    global_history.row(0) = C.transpose();

    for (int k = 1; k <= steps; ++k) {
        double tn = t;
        double tn1 = t + h;
        Eigen::VectorXd fn = f(tn, C);

    switch (integrator_) {
            case Integrator::rectangle_gauche: {
                C = C + Integration::rectangle_gauche(fn, h);
                break;
            }
            case Integrator::rectangle_droite: {
                // predictor value f at tn1 using explicit Euler predictor
                Eigen::VectorXd Cp = C + Integration::rectangle_gauche(fn, h);
                Eigen::VectorXd fn1 = f(tn1, Cp);
                C = C + Integration::rectangle_droite(fn1, h);
                break;
            }
            case Integrator::trapeze: {
                // implicit trapezoidal via fixed-point on G(C1) = C0 + 0.5*h*(f(tn,C0) + f(tn1,C1))
                Eigen::VectorXd C0 = C;
                auto G = [&](const Eigen::VectorXd& C1)->Eigen::VectorXd {
                    return C0 + 0.5 * h * (fn + f(tn1, C1));
                };
                t += h;

                // After computing new state C, update compartment quantities and append history
                if (C.size() >= 4) {
                    atmosphere.setQuantite(C(0));
                    arbre.setQuantite(C(1));
                    sol.setQuantite(C(2));
                    ocean.setQuantite(C(3));

                    atmosphere.appendHistory(C(0));
                    arbre.appendHistory(C(1));
                    sol.appendHistory(C(2));
                    ocean.appendHistory(C(3));
                }

                global_history.row(k) = C.transpose();
            case Integrator::pointFixe: {
                // user-supplied fixed-point iteration on map F(C) = C + h*f(tn,C)
                auto F = [&](const Eigen::VectorXd& X)->Eigen::VectorXd {
                    return X + h * f(tn, X);
                };
                C = Resolution::pointFixe(C, F);
                break;
            }
            case Integrator::newton: {
                // use Newton to solve G(C1) = C1 - C0 - h*f(tn1, C1) = 0 (backward Euler)
                Eigen::VectorXd C0 = C;
                auto G = [&](const Eigen::VectorXd& X)->Eigen::VectorXd {
                    return X - C0 - h * f(tn1, X);
                };
                // approximate Jacobian via finite differences
                auto dG = [&](const Eigen::VectorXd& X)->Eigen::MatrixXd {
                    const double eps_fd = 1e-8;
                    int n = X.size();
                    Eigen::MatrixXd J(n,n);
                    Eigen::VectorXd FX = f(tn1, X);
                    for (int j = 0; j < n; ++j) {
                        Eigen::VectorXd XP = X;
                        XP(j) += eps_fd;
                        Eigen::VectorXd FXP = f(tn1, XP);
                        Eigen::VectorXd diff = (FXP - FX) / eps_fd;
                        J.col(j) = -h * diff;
                    }
                    // add identity for d(X - C0)/dX
                    J += Eigen::MatrixXd::Identity(n,n);
                    return J;
                };
                C = Resolution::newton(C, G, dG);
                break;
            }
        }
        t += h;
        global_history.row(k) = C.transpose();
    }
    current_pos = steps;
}

void Simulation::exportHistoryCSV(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("could not open file for writing: " + filename);

    // header: time, state0, state1, ...
    out << "time";
    for (int j = 0; j < global_history.cols(); ++j) out << ",C" << j;
    out << '\n';

    for (int i = 0; i < global_history.rows(); ++i) {
        double time = t0 + i * h;
        out << std::fixed << std::setprecision(8) << time;
        for (int j = 0; j < global_history.cols(); ++j) {
            out << ',' << std::scientific << std::setprecision(12) << global_history(i,j);
        }
        out << '\n';
    }
}

int Simulation::getMaxIter() const { return max_iter; }
void Simulation::setMaxIter(int max_iter_) { max_iter = max_iter_; }

double Simulation::getTolerance() const { return eps; }
void Simulation::setTolerance(double eps_) { eps = eps_; }

double Simulation::getStepSize() const { return h; }
void Simulation::setStepSize(double h_) { if (h_ > 0) h = h_; }

const Eigen::VectorXd& Simulation::getCondInitial() const { return C0; }
void Simulation::setCondInitial(const Eigen::VectorXd& C0_) { C0 = C0_; C = C0_; }

double Simulation::getSimulBeginTime() const { return t0; }
void Simulation::setSimulBeginTime(double t0_) { t0 = t0_; }

double Simulation::getSimulTime() const { return t0 + current_pos * h; }
void Simulation::setSimulTime(double t) { t0 = t; }

double Simulation::getSimulEndTime() const { return tf; }
void Simulation::setSimulEndTime(double tf_) { tf = tf_; }

const Eigen::MatrixXd& Simulation::getSimulHistory() const { return global_history; }

const Eigen::VectorXd& Simulation::getSystemState() const { return C; }
void Simulation::setSystemState(const Eigen::VectorXd& C_) { C = C_; }