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
#include "csv.hpp"

void plot();

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    constexpr double pi = std::numbers::pi;
    std::vector<int> ns = {8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64};
    std::vector<double> serial_times, omp_times, mpi_times, hybrid_times;
    std::vector<double> omp_speedups, mpi_speedups, hybrid_speedups;
    std::vector<double> l2_errors;

    // Only print headers on rank 0
    if (rank == 0)
    {
        std::cout << std::setw(8) << "n"
                  << std::setw(15) << "Serial Time(s)"
                  << std::setw(15) << "OMP Time(s)"
                  << std::setw(15) << "MPI Time(s)"
                  << std::setw(15) << "Hybrid Time(s)"
                  << std::setw(10) << "OMP SU"
                  << std::setw(10) << "MPI SU"
                  << std::setw(10) << "Hybrid SU"
                  << std::setw(15) << "L2 error" << "\n";
        std::cout << "----------------------------------------------------------------------------------------------------------------------\n";
    }

    for (int n : ns)
    {
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

        double serial_time = 0.0, omp_time = 0.0, mpi_time = 0.0, hybrid_time = 0.0;
        double serial_l2 = 0.0;

        // Only run serial and OpenMP tests on rank 0 to avoid duplication
        if (rank == 0)
        {
            // Serial test
            auto start = std::chrono::high_resolution_clock::now();
            solver.solve_serial();
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> serial_elapsed = end - start;
            serial_time = serial_elapsed.count();
            serial_l2 = solver.l2_error();

            solver.reset(); // Reset the solver for the next run

            // OpenMP test
            start = std::chrono::high_resolution_clock::now();
            solver.solve_omp();
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> omp_elapsed = end - start;
            omp_time = omp_elapsed.count();
        }

        // Reset solver for MPI run
        solver.reset();

        // MPI test (all processes participate)
        auto start_mpi = std::chrono::high_resolution_clock::now();
        solver.solve_mpi();
        auto end_mpi = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> mpi_elapsed = end_mpi - start_mpi;
        mpi_time = mpi_elapsed.count();

        // Reset solver for hybrid run
        solver.reset();

        // Hybrid OpenMP+MPI test (all processes participate)
        auto start_hybrid = std::chrono::high_resolution_clock::now();
        solver.solve_hybrid();
        auto end_hybrid = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> hybrid_elapsed = end_hybrid - start_hybrid;
        hybrid_time = hybrid_elapsed.count();

        // Only rank 0 handles output and data collection
        if (rank == 0)
        {
            double omp_speedup = serial_time / omp_time;
            double mpi_speedup = serial_time / mpi_time;
            double hybrid_speedup = serial_time / hybrid_time;

            serial_times.push_back(serial_time);
            omp_times.push_back(omp_time);
            mpi_times.push_back(mpi_time);
            hybrid_times.push_back(hybrid_time);
            omp_speedups.push_back(omp_speedup);
            mpi_speedups.push_back(mpi_speedup);
            hybrid_speedups.push_back(hybrid_speedup);
            l2_errors.push_back(serial_l2);

            std::cout << std::setw(8) << n
                      << std::setw(15) << std::fixed << std::setprecision(6) << serial_time
                      << std::setw(15) << std::fixed << std::setprecision(6) << omp_time
                      << std::setw(15) << std::fixed << std::setprecision(6) << mpi_time
                      << std::setw(15) << std::fixed << std::setprecision(6) << hybrid_time
                      << std::setw(10) << std::fixed << std::setprecision(4) << omp_speedup
                      << std::setw(10) << std::fixed << std::setprecision(4) << mpi_speedup
                      << std::setw(10) << std::fixed << std::setprecision(4) << hybrid_speedup
                      << std::setw(15) << std::scientific << std::setprecision(3) << serial_l2 << "\n";

            if (n == 64)
                solver.save_vtk("solution_" + std::to_string(size) + "_n_" + std::to_string(n));
        }
    }

    // Only write results file on rank 0
    if (rank == 0)
    {
        std::ofstream ofs("test/data/results_" + std::to_string(size) + ".csv");
        ofs << "n,serial,omp,mpi,hybrid,omp_speedup,mpi_speedup,hybrid_speedup,l2_error\n";
        for (size_t i = 0; i < ns.size(); ++i)
        {
            ofs << ns[i] << "," << serial_times[i] << "," << omp_times[i] << "," << mpi_times[i] << "," << hybrid_times[i]
                << "," << omp_speedups[i] << "," << mpi_speedups[i] << "," << hybrid_speedups[i] << "," << l2_errors[i] << "\n";
        }
        ofs.close();
    }

    if (size == 8 && rank == 0)
    {
        plot();
    }

    MPI_Finalize();
    return 0;
}

void plot()
{
    try
    {
        scalabilityTest();
        gridSizeTest("results_1.csv");
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return;
}