#include "simulation/simulation.h"
#include <Eigen/Dense>
#include <iostream>

// simple RHS: dC/dt = -0.1 * C
Eigen::VectorXd rhs(double t, const Eigen::VectorXd& C) {
    return -0.1 * C;
}

int main() {
    Eigen::VectorXd C0(3);
    C0 << 100.0, 50.0, 10.0;

    Simulation sim(C0, 0.0, 10.0, 1.0);
    sim.setSimulationMethod(rhs);
    sim.run();
    std::cout << "Final state:\n" << sim.getSystemState() << std::endl;
    sim.exportHistoryCSV("sim_history.csv");
    return 0;
}
