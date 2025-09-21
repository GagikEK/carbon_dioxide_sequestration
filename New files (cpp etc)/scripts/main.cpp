#include "simulation/simulation.h"
#include <Eigen/Dense>
#include <iostream>

int main() {
     // --- Simulation ---
        Simulation simul;
        simul.setStepSize(0.01);
        simul.setSolver(Simulation::Solver::EulerExplicite);
        simul.run();

        // --- Récupérer et afficher l'état final ---

        simul.exportHistoryTXT("simulation_history2.txt");
        simul.exportHistoryCSV("simulation_history2.csv");

    return 0;
}
