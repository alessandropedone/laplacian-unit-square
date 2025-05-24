#include <iostream>
#include <vector>
#include <cmath>
#include <numbers>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <omp.h>
#include <mpi.h>

#include "jacobi_solver.hpp"
#include "vtk.hpp"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    constexpr double pi = std::numbers::pi;
    std::vector<int> ns = {32}; //{8, 16, 32, 64, 128};
    std::vector<double> serial_times, omp_times, mpi_times, omp_speedups, mpi_speedups;
    std::vector<double> l2_errors; // Fixed typo: was "l2_errrors"

    // Only print headers on rank 0
    if (rank == 0) {
        std::cout << std::setw(8) << "n"
                  << std::setw(20) << "Serial Time (s)"
                  << std::setw(20) << "OMP Time (s)"
                  << std::setw(20) << "MPI Time (s)"
                  << std::setw(15) << "OMP Speedup"
                  << std::setw(15) << "MPI Speedup"
                  << std::setw(20) << "L2 error" << "\n";
        std::cout << "---------------------------------------------------------------------------------------------------------------------\n";
    }

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

        double serial_time = 0.0, omp_time = 0.0, mpi_time = 0.0;
        double serial_l2 = 0.0;

        // Only run serial and OpenMP tests on rank 0 to avoid duplication
        if (rank == 0) {
            auto start = std::chrono::high_resolution_clock::now();
            solver.solve_serial();
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> serial_elapsed = end - start;
            serial_time = serial_elapsed.count();
            serial_l2 = solver.l2_error();

            solver.reset(); // Reset the solver for the next run

            start = std::chrono::high_resolution_clock::now();
            solver.solve_omp();
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> omp_elapsed = end - start;
            omp_time = omp_elapsed.count();

            if (n == 128) {
                solver.save_vtk("solution");
            }
        }

        // Reset solver for MPI run
        solver.reset();
        
        // Run MPI solver (all processes participate)
        auto start_mpi = std::chrono::high_resolution_clock::now();
        solver.solve_mpi();
        auto end_mpi = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> mpi_elapsed = end_mpi - start_mpi;
        mpi_time = mpi_elapsed.count();

        // Only rank 0 handles output and data collection
        if (rank == 0) {
            double omp_speedup = serial_time / omp_time;
            double mpi_speedup = serial_time / mpi_time;

            serial_times.push_back(serial_time);
            omp_times.push_back(omp_time);
            mpi_times.push_back(mpi_time);
            omp_speedups.push_back(omp_speedup);
            mpi_speedups.push_back(mpi_speedup);
            l2_errors.push_back(serial_l2);

            std::cout << std::setw(8) << n
                      << std::setw(20) << serial_time
                      << std::setw(20) << omp_time
                      << std::setw(20) << mpi_time
                      << std::setw(15) << omp_speedup
                      << std::setw(15) << mpi_speedup
                      << std::setw(20) << serial_l2 << "\n";
        }
    }

    // Only write results file on rank 0
    if (rank == 0) {
        std::ofstream ofs("results.csv");
        ofs << "n,serial,omp,mpi,omp_speedup,mpi_speedup,l2_error\n";
        for (size_t i = 0; i < ns.size(); ++i) {
            ofs << ns[i] << "," << serial_times[i] << "," << omp_times[i] << "," << mpi_times[i]
                << "," << omp_speedups[i] << "," << mpi_speedups[i] << "," << l2_errors[i] << "\n";
        }
        ofs.close();
    }

    MPI_Finalize();
    return 0;
}