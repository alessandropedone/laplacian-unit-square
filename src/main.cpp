
#include <iostream>
#include <vector>
#include <cmath>
#include <numbers>
#include <iomanip> // for std::setw

#include "jacobi_serial_solver.hpp"
#include "vtk.hpp"



int main() {
    constexpr double pi = std::numbers::pi;

    int n = 100;

    JacobiSerialSolver solver (
        std::vector<double>(n * n, 0.0), // initial guess
        [=](double x, double y) { return 8 * pi * pi * sin(2 * pi * x) * sin(2 * pi * y); }, // rhs
        [=](double x, double y) { return 0.0; }, // top boundary condition
        [=](double x, double y) { return 0.0; }, // right boundary condition
        [=](double x, double y) { return 0.0; }, // bottom boundary condition
        [=](double x, double y) { return 0.0; }, // left boundary condition
        n, // grid size
        40000, // max iterations
        1e-15, // tolerance
        [=](double x, double y) { return sin(2 * pi * x) * sin(2 * pi * y); } // exact solution
    );

    
    
    solver.solve();

    return 0;
}