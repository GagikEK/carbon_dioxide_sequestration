#include <Eigen/Dense>
#include <iostream>

int main() {
    // Définir une matrice 2x2
    Eigen::Matrix2d mat;
    mat << 1, 2,
           3, 4;

    // Définir un vecteur 2x1
    Eigen::Vector2d vec(5, 6);

    // Multiplication matrice-vecteur
    Eigen::Vector2d result = mat * vec;

    std::cout << "Résultat de mat * vec :\n" << result << std::endl;

    return 0;
}
