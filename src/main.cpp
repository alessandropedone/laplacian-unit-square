
#include <iostream>
#include <vector>
#include <cmath>
#include <numbers>

#include "serial_solver.hpp"
#include "vtk.hpp"

int main() {
    constexpr double pi = std::numbers::pi;

    const size_t n = 5;
    /// @brief grid points over [0, 1]
    std::vector<double> x_points(n), y_points(n);
    for (size_t i = 0; i < n; ++i) {
        x_points[i] = i * 1.0 / n;
        y_points[i] = x_points[i];
    }

    std::vector<double> exact_sol(n * n, 0.0), top_bc(n, 0.0), bottom_bc(n, 0.0), right_bc(n-2, 0.0), left_bc(n-2, 0.0);

    // assign_exact_sol
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            exact_sol[i * n + j] = sin(2 * pi * x_points[i]) * sin(2 * pi * y_points[j]);
        }
    }

    std::vector<double> initial_guess(n*n, 0.0), rhs(n*n, 0.0);

    // assign_rhs
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            rhs[i * n + j] = 8 * pi * pi * sin(2 * pi * x_points[i]) * sin(2 * pi * y_points[j]);
        }
    }

    unsigned max_iter = 1000;
    double tol = 1e-6;

    SerialSolver solver(exact_sol, initial_guess, rhs, top_bc, right_bc, bottom_bc, left_bc, n, max_iter, tol);

    solver.solve(x_points, y_points);
    return 0;
}