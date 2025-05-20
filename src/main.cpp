
#include <iostream>
#include <vector>
#include <cmath>
#include <numbers>
#include <iomanip>
#include <omp.h>
#include <mpi.h>

#include "jacobi_serial_solver.hpp"
#include "vtk.hpp"

int main(int argc, char *argv[])
{
    constexpr double pi = std::numbers::pi;

    int n = 100;

    JacobiSerialSolver solver(
        std::vector<double>(n * n, 0.0), // initial guess
        [=](double x, double y)
        { return 8 * pi * pi * sin(2 * pi * x) * sin(2 * pi * y); }, // rhs
        [=](double x, double y)
        { return 0.0; }, // top boundary condition
        [=](double x, double y)
        { return 0.0; }, // right boundary condition
        [=](double x, double y)
        { return 0.0; }, // bottom boundary condition
        [=](double x, double y)
        { return 0.0; }, // left boundary condition
        n,               // grid size
        10000,              // max iterations
        1e-10,           // tolerance
        [=](double x, double y)
        { return sin(2 * pi * x) * sin(2 * pi * y); } // exact solution
    );

    solver.solve();

    return 0;
}