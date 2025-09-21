cmake_minimum_required(VERSION 3.16)
project(SimulationCarbone LANGUAGES CXX)

# Utiliser C++17 (ou C++20 si besoin)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Trouver Eigen (assume qu'il est installé via apt / brew ou trouvé par CMake)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Liste tous tes fichiers sources
set(SOURCES
    main.cpp
    simulation/simulation.cpp
    compartiments/Compartiments.cpp
    compartiments/Arbres.cpp
    compartiments/Sol.cpp
    compartiments/Atmosphere.cpp
    compartiments/Humains.cpp
    compartiments/Oceans.cpp
    methodes/integration.cpp
    methodes/resolution.cpp
)

# Crée l'exécutable
add_executable(main ${SOURCES})

# Lier Eigen
target_link_libraries(main PRIVATE Eigen3::Eigen)