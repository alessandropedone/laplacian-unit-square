
#include <iostream>
#include <vector>
#include <cmath>
#include <numbers>
#include <iomanip>
#include <chrono>
#include <omp.h>
#include <mpi.h>

#include "jacobi_solver.hpp"
#include "vtk.hpp"

int main(int argc, char *argv[])
{
    constexpr double pi = std::numbers::pi;
    std::vector<int> ns = {8, 16, 32, 64, 128};
    std::vector<double> serial_times, omp_times, speedups;
    std::vector<double> l2_errrors;

    std::cout << std::setw(8) << "n"
              << std::setw(20) << "Serial Time (s)"
              << std::setw(20) << "OMP Time (s)"
              << std::setw(15) << "Speedup"
              << std::setw(20) << "L2 error"<< "\n";
    std::cout << "---------------------------------------------------------------------------------------------\n";

    for (int n : ns) {
        JacobiSolver solver(
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
            30000,           // max iterations
            1e-15,           // tolerance
            [=](double x, double y)
            { return sin(2 * pi * x) * sin(2 * pi * y); } // exact solution
        );

        auto start = std::chrono::high_resolution_clock::now();
        solver.solve_serial();
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> serial_elapsed = end - start;
        double serial_l2 = solver.l2_error();

        solver.reset(); // Reset the solver for the next run

        start = std::chrono::high_resolution_clock::now();
        solver.solve_omp();
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> omp_elapsed = end - start;

        double speedup = serial_elapsed.count() / omp_elapsed.count();

        serial_times.push_back(serial_elapsed.count());
        omp_times.push_back(omp_elapsed.count());
        speedups.push_back(speedup);
        l2_errrors.push_back(serial_l2);

        std::cout << std::setw(8) << n
                  << std::setw(20) << serial_elapsed.count()
                  << std::setw(20) << omp_elapsed.count()
                  << std::setw(15) << speedup
                  << std::setw(20) << serial_l2 << "\n";

        if (n == 128) {
            solver.save_vtk("solution");
        }
        
    }

    // Output data for plotting
    std::ofstream ofs("results.csv");
    ofs << "n,serial,omp,speedup,l2_error\n";
    for (size_t i = 0; i < ns.size(); ++i) {
        ofs << ns[i] << "," << serial_times[i] << "," << omp_times[i] << "," << speedups[i]
            << "," << l2_errrors[i] << "\n";
    }
    ofs.close();
    return 0;
}