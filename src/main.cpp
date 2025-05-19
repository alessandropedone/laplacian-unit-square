
#include <iostream>
#include <vector>
#include <cmath>
#include <numbers>
#include <iomanip> // for std::setw

#include "serial_solver.hpp"
#include "vtk.hpp"

int main() {
    constexpr double pi = std::numbers::pi;

    const size_t n = 10;
    /// @brief grid points over [0, 1]
    std::vector<double> x_points(n), y_points(n);
    for (size_t i = 0; i < n; ++i) {
        x_points[i] = i * 1.0 / (n-1);
        y_points[i] = x_points[i];
    }

    std::vector<double> exact_sol(n * n, 0.0), topbc(n, 0.0), bottombc(n, 0.0), rightbc(n, 0.0), leftbc(n, 0.0);

    // assign exact_sol
    std::cout << "Exact solution:" << std::endl;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            exact_sol[i * n + j] = sin(2 * pi * x_points[i]) * sin(2 * pi * y_points[j]);
            std::cout << std::setw(10) << exact_sol[i * n + j] << " ";
        }
        std::cout << std::endl;
    }

    std::vector<double> initial_guess(n*n, 0.), rhs(n*n, 0.0);

    // assign rhs
    std::cout << "RHS:" << std::endl;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            rhs[i * n + j] = 8 * pi * pi * sin(2 * pi * x_points[i]) * sin(2 * pi * y_points[j]);
            std::cout << std::setw(10) << rhs[i * n + j] << " ";
        }
        std::cout << std::endl;
    }

    unsigned max_iter = 100000;
    double tol = 1e-18;

    SerialSolver solver(exact_sol, initial_guess, rhs, topbc, rightbc, bottombc, leftbc, n, max_iter, tol);

    solver.solve(x_points, y_points);
    return 0;
}