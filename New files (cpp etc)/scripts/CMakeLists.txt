cmake_minimum_required(VERSION 3.10)
project(Exercice3 LANGUAGES CXX)

# Use C++17 (Eigen works with older standards, but C++17 is recommended)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Try to find an installed Eigen (preferred)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(mon_exe main.cpp)

# Link Eigen (header-only target sets include directories)
target_link_libraries(mon_exe PRIVATE Eigen3::Eigen)


