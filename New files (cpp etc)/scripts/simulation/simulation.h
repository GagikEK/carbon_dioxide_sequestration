#ifndef Simulation_H
#define Simulation_H

#include <Eigen/Dense>
#include <functional>
#include <vector>

// compartment headers (relative path from scripts/simulation)
#include "../compartiments/Compartiments.h"

// base/default right-hand-side function signature
Eigen::VectorXd f_base(double tn, const Eigen::VectorXd& Cn);

class Simulation {
private:
    // Concrete compartment objects (use the plural class names defined in headers)
    Arbres arbre;
    Oceans ocean;
    Humains humain;
    Sol sol;
    Atmosphere atmosphere;

    // Numerical / control parameters
    int max_iter{10000};        // maximum iterations for nonlinear solves (if used)
    double eps{1e-6};           // tolerance for iterative solvers
    double h{1.0};              // timestep size

    // time / initial conditions
    Eigen::VectorXd C0;         // initial state vector
    double t0{2024.0};
    double tf{2100.0};

    // history storage: rows = time steps, cols = state size
    Eigen::MatrixXd global_history;
    Eigen::VectorXd C;         // current system state
    int current_pos{0};
    Eigen::VectorXd time_vector; // store times corresponding to history rows

    // user-provided right-hand-side function (takes time and state by const-ref)
    std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)> f;

    // helper: number of time steps (computed from t0, tf and h)
    int nb_pas() const;

public:
    enum class Solver { EulerExplicite, RectangleDroite_PointFixe, TrapezePointFixe, Trapeze_Newton, RectangleDroite_Newton };

private:
    // selected Solver (declared after enum)
    Solver Solver_{Solver::EulerExplicite};

public:
    // ctors / dtor
    Simulation();
    Simulation(const Eigen::VectorXd& C0,
               double t0,
               double tf,
               double h,
               int max_iter = 10000,
               double eps = 1e-6,
               const Arbres& arbre = Arbres(),
               const Sol& sol = Sol(),
               const Atmosphere& atmosphere = Atmosphere(),
               const Humains& humain = Humains(),
               const Oceans& ocean = Oceans());
    ~Simulation();

    // run the simulation (will fill global_history)
    void run();

    // Export history to CSV
    void exportHistoryCSV(const std::string& filename) const;

    void exportHistoryTXT(const std::string& filename) const; 

    // getters / setters (const-correct, avoid copies where sensible)
    int getMaxIter() const;
    void setMaxIter(int max_iter);

    double getTolerance() const;
    void setTolerance(double eps);

    double getStepSize() const;
    void setStepSize(double h);

    const Eigen::VectorXd& getCondInitial() const;
    void setCondInitial(const Eigen::VectorXd& C0);

    double getSimulBeginTime() const;
    void setSimulBeginTime(double t0);

    double getSimulTime() const; // current simulation time
    void setSimulTime(double t);

    double getSimulEndTime() const;
    void setSimulEndTime(double tf);

    const Eigen::MatrixXd& getSimulHistory() const;

    const Eigen::VectorXd& getSystemState() const;
    void setSystemState(const Eigen::VectorXd& C);

    // convenience: set RHS using raw function pointer too
    void setSimulFunction(Eigen::VectorXd (*fp)(double, const Eigen::VectorXd&));

    // set/get Solver
    void setSolver(Solver Solver_);
    Solver getSolver() const;
};

#endif