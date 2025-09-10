#include <iostream>
#include <vector>
#include <fstream>
#include "compartiments.h"
#include "Atmosphere.h"
#include "Arbres.h"
#include "Sol.h"
#include "Oceans.h"
#include "Humains.h"

int main() {
    std::cout << "=== Carbon Dioxide Sequestration Simulation ===" << std::endl;
    
    // Simulation parameters
    const int simulation_years = 100;
    const int history_size = simulation_years + 1;
    
    // Initialize compartments with default values
    // You can modify these values based on your specific simulation needs
    double alpha = 0.04;
    double beta = 0.02;
    double gamma = 0.015;
    double delta = 0.008;
    double K = 900;
    
    // Set global constants
    Compartiments::setConstantes(alpha, beta, gamma, delta, K);
    
    // Initialize compartments
    Atmosphere atmosphere(800, 0, alpha, beta, gamma, delta, K, history_size);
    Arbres arbres(500, 0, alpha, beta, gamma, delta, K, history_size);
    Sol sol(2000, 0, alpha, beta, gamma, delta, K, history_size);
    Oceans oceans(38000, 0, alpha, beta, gamma, delta, K, history_size);
    Humains humains(0, 0, alpha, beta, gamma, delta, K, history_size, 0.1);
    
    std::cout << "Initial values:" << std::endl;
    std::cout << "Atmosphere: " << atmosphere.getQuantite() << " GtC" << std::endl;
    std::cout << "Trees: " << arbres.getQuantite() << " GtC" << std::endl;
    std::cout << "Soil: " << sol.getQuantite() << " GtC" << std::endl;
    std::cout << "Oceans: " << oceans.getQuantite() << " GtC" << std::endl;
    std::cout << "Humans: " << humains.getQuantite() << " GtC" << std::endl;
    std::cout << std::endl;
    
    // Store initial values in history
    atmosphere.appendHistory(atmosphere.getQuantite());
    arbres.appendHistory(arbres.getQuantite());
    sol.appendHistory(sol.getQuantite());
    oceans.appendHistory(oceans.getQuantite());
    humains.appendHistory(humains.getQuantite());
    
    // Simulation loop
    std::cout << "Starting simulation for " << simulation_years << " years..." << std::endl;
    
    for (int year = 1; year <= simulation_years; year++) {
        // Get current values
        double CA = atmosphere.getQuantite();
        double CT = arbres.getQuantite();
        double CS = sol.getQuantite();
        double CO = oceans.getQuantite();
        double CH = humains.getQuantite();
        
        // Update each compartment
        double new_CA = atmosphere.update(CA, CT, CS);
        double new_CT = arbres.update(CA, CT, CS);
        double new_CS = sol.update(arbres, sol, atmosphere, humains, oceans);
        double new_CO = oceans.update(CA, CT, CS);
        double new_CH = humains.update(arbres, sol, atmosphere, humains, oceans);
        
        // Update quantities
        atmosphere.setQuantite(new_CA);
        arbres.setQuantite(new_CT);
        sol.setQuantite(new_CS);
        oceans.setQuantite(new_CO);
        humains.setQuantite(new_CH);
        
        // Store in history
        atmosphere.appendHistory(new_CA);
        arbres.appendHistory(new_CT);
        sol.appendHistory(new_CS);
        oceans.appendHistory(new_CO);
        humains.appendHistory(new_CH);
        
        // Print progress every 10 years
        if (year % 10 == 0) {
            std::cout << "Year " << year << ":" << std::endl;
            std::cout << "  Atmosphere: " << new_CA << " GtC" << std::endl;
            std::cout << "  Trees: " << new_CT << " GtC" << std::endl;
            std::cout << "  Soil: " << new_CS << " GtC" << std::endl;
            std::cout << "  Oceans: " << new_CO << " GtC" << std::endl;
            std::cout << "  Humans: " << new_CH << " GtC" << std::endl;
            std::cout << std::endl;
        }
    }
    
    std::cout << "Simulation completed!" << std::endl;
    std::cout << "Final values:" << std::endl;
    std::cout << "Atmosphere: " << atmosphere.getQuantite() << " GtC" << std::endl;
    std::cout << "Trees: " << arbres.getQuantite() << " GtC" << std::endl;
    std::cout << "Soil: " << sol.getQuantite() << " GtC" << std::endl;
    std::cout << "Oceans: " << oceans.getQuantite() << " GtC" << std::endl;
    std::cout << "Humans: " << humains.getQuantite() << " GtC" << std::endl;
    
    return 0;
}
