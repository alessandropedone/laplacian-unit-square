/**
 * @file main.cpp
 * @brief Main driver program for parallel Jacobi solver performance testing and analysis.
 *
 * This program implements a comprehensive performance analysis of different parallelization
 * strategies for solving the 2D Poisson equation using the Jacobi iterative method. It tests
 * serial, OpenMP, MPI, hybrid (OpenMP+MPI), and direct solver implementations across various 
 * grid sizes.
 *
 * The program solves the equation:
 * -∇²u = 8π²sin(2πx)sin(2πy)
 * with homogeneous Dirichlet boundary conditions on the unit square [0,1]×[0,1].
 * The exact solution is u(x,y) = sin(2πx)sin(2πy).
 *
 * Features:
 * - Performance benchmarking across multiple grid sizes (8×8 to 64×64)
 * - Five different solver implementations: serial, OpenMP, MPI, hybrid, and direct MPI
 * - Speedup calculations for each parallelization method
 * - L2 error analysis for accuracy validation
 * - VTK output for visualization of the largest grid solution
 * - CSV data export for further analysis
 * - Automated plotting and scalability analysis
 * - Optional parameter file support via GetPot and muParserX
 *
 * Command Line Options:
 * - --use-datafile or -d: Read parameters from data.txt file (slower due to parser overhead)
 *
 * Output:
 * - Console table showing execution times, speedups, and errors for all methods
 * - CSV files with detailed results for each MPI process count
 * - VTK files for solution visualization (n=64)
 * - Performance plots (when run with 4 MPI processes)
 *
 * @note This program requires MPI initialization and should be run with multiple processes
 *       to evaluate MPI and hybrid performance. Only rank 0 handles output operations.
 *
 * @note Parameter file reading is available but introduces overhead due to muParserX interface.
 */
#include <iostream>
#include <vector>
#include <cmath>
#include <numbers>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <omp.h>
#include <mpi.h>
#include <GetPot>

#include "muparser_interface.hpp"
#include "solver.hpp"
#include "vtk.hpp"
#include "plot.hpp"
#include "simulation_parameters.hpp"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Possibility to read the parameters from a file but the test runs a lot slower
    // because of the overhead of muparserx interface.
    bool use_datafile = false;
    if (argc > 1)
    {
        std::string arg = argv[1];
        if (arg == "--use-datafile" || arg == "-d")
        {
            use_datafile = true;
        }
    }

    solver::SimulationParameters params;

    if (use_datafile)
    {
        // Read parameters from data.txt file (only on root process)
        if (rank == 0)
        {
            const GetPot datafile("data.txt");
            params.f_str = datafile("f", "8 * pi * pi * sin(2 * pi * x[0]) * sin(2 * pi * x[1])");
            params.uex_str = datafile("uex", "sin(2 * pi * x[0]) * sin(2 * pi * x[1])");
            params.bc_top_str = datafile("d_bc_top", "0.0");
            params.bc_right_str = datafile("d_bc_right", "0.0");
            params.bc_bottom_str = datafile("d_bc_bottom", "0.0");
            params.bc_left_str = datafile("d_bc_left", "0.0");
            params.tol = datafile("tol", 1e-15);
            params.max_iter = datafile("max_iter", 30000);
        }

        // Broadcast all parameters to all processes
        params.broadcast(0, MPI_COMM_WORLD);
    }

    std::vector<int> ns = {8, 16, 24, 32, 40, 48, 56, 64};
    std::vector<double> serial_times, omp_times, mpi_times, hybrid_times, direct_times;
    std::vector<double> omp_speedups, mpi_speedups, hybrid_speedups, direct_speedups;
    std::vector<double> l2_errors;

    // Only print headers on rank 0
    if (rank == 0)
    {
        std::cout << std::setw(8) << "n"
                  << std::setw(15) << "Serial Time(s)"
                  << std::setw(15) << "OMP Time(s)"
                  << std::setw(15) << "MPI Time(s)"
                  << std::setw(15) << "Hybrid Time(s)"
                  << std::setw(15) << "Direct Time(s)"
                  << std::setw(10) << "OMP SU"
                  << std::setw(10) << "MPI SU"
                  << std::setw(10) << "Hybrid SU"
                  << std::setw(10) << "Direct SU"
                  << std::setw(15) << "L2 error" << "\n";
        std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------\n";
    }

    for (int n : ns)
    {

        // Possibility to read the parameters from a file but the test runs a lot slower
        // because of the overhead of muparserx interface.
        solver::Solver solver;
        constexpr auto pi = std::numbers::pi;
        if (use_datafile)
        {
            // Create muParserX interfaces
            muparser::muParserXInterface f(params.f_str, 2);
            muparser::muParserXInterface uex(params.uex_str, 2);
            muparser::muParserXInterface top_bc(params.bc_top_str, 2);
            muparser::muParserXInterface right_bc(params.bc_right_str, 2);
            muparser::muParserXInterface bottom_bc(params.bc_bottom_str, 2);
            muparser::muParserXInterface left_bc(params.bc_left_str, 2);
            solver.set_bc(top_bc, right_bc, bottom_bc, left_bc);       // Set boundary conditions
            solver.set_initial_guess(std::vector<double>(n * n, 0.0)); // Initial guess
            solver.set_f(f);                                           // Set right-hand side function
            solver.set_uex(uex);                                       // Set exact solution function
            solver.set_n(n);                                           // Set grid size
            solver.set_max_iter(params.max_iter);                      // Set maximum iterations
            solver.set_tol(params.tol);                                // Set tolerance for convergence
        }
        else
        {

            solver.set_bc([=](std::vector<double> x)
                          { return 0.0; }, [=](std::vector<double> x)
                          { return 0.0; }, [=](std::vector<double> x)
                          { return 0.0; }, [=](std::vector<double> x)
                          { return 0.0; });                            // Set boundary conditions
            solver.set_initial_guess(std::vector<double>(n * n, 0.0)); // Initial guess
            solver.set_f([=](std::vector<double> x)
                         { return 8 * pi * pi * sin(2 * pi * x[0]) * sin(2 * pi * x[1]); }); // Set right-hand side function
            solver.set_uex([=](std::vector<double> x)
                           { return sin(2 * pi * x[0]) * sin(2 * pi * x[1]); }); // Set exact solution function
            solver.set_n(n);                                                     // Set grid size
            solver.set_max_iter(30000);                                          // Set maximum iterations
            solver.set_tol(1e-15);                                               // Set tolerance for convergence
        }

        double serial_time = 0.0, omp_time = 0.0, mpi_time = 0.0, hybrid_time = 0.0, direct_time = 0.0;
        double serial_l2 = 0.0;

        // Only run serial and OpenMP tests on rank 0 to avoid duplication
        if (rank == 0)
        {
            // Serial test
            auto start = std::chrono::high_resolution_clock::now();
            solver.solve_jacobi_serial();
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> serial_elapsed = end - start;
            serial_time = serial_elapsed.count();
            serial_l2 = solver.l2_error();

            // Reset the solver for OMP run
            solver.reset();

            // OpenMP test
            start = std::chrono::high_resolution_clock::now();
            solver.solve_jacobi_omp();
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> omp_elapsed = end - start;
            omp_time = omp_elapsed.count();
        }

        // Reset solver for MPI run
        solver.reset();

        // MPI test (all processes participate)
        auto start_mpi = std::chrono::high_resolution_clock::now();
        solver.solve_jacobi_mpi();
        auto end_mpi = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> mpi_elapsed = end_mpi - start_mpi;
        mpi_time = mpi_elapsed.count();

        // Reset solver for hybrid run
        solver.reset();

        // Hybrid OpenMP+MPI test (all processes participate)
        auto start_hybrid = std::chrono::high_resolution_clock::now();
        solver.solve_jacobi_hybrid();
        auto end_hybrid = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> hybrid_elapsed = end_hybrid - start_hybrid;
        hybrid_time = hybrid_elapsed.count();

        // Reset solver for direct local solver test
        solver.reset();

        // Mpi test with direct local solver
        auto direct_start = std::chrono::high_resolution_clock::now();
        solver.solve_direct_mpi();
        auto direct_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> direct_elapsed = direct_end - direct_start;
        direct_time = direct_elapsed.count();

        // Only rank 0 handles output and data collection
        if (rank == 0)
        {
            double omp_speedup = serial_time / omp_time;
            double mpi_speedup = serial_time / mpi_time;
            double hybrid_speedup = serial_time / hybrid_time;
            double direct_speedup = serial_time / direct_time;

            serial_times.push_back(serial_time);
            omp_times.push_back(omp_time);
            mpi_times.push_back(mpi_time);
            hybrid_times.push_back(hybrid_time);
            direct_times.push_back(direct_time);
            omp_speedups.push_back(omp_speedup);
            mpi_speedups.push_back(mpi_speedup);
            hybrid_speedups.push_back(hybrid_speedup);
            direct_speedups.push_back(direct_speedup);
            l2_errors.push_back(serial_l2);

            std::cout << std::setw(8) << n
                      << std::setw(15) << std::fixed << std::setprecision(6) << serial_time
                      << std::setw(15) << std::fixed << std::setprecision(6) << omp_time
                      << std::setw(15) << std::fixed << std::setprecision(6) << mpi_time
                      << std::setw(15) << std::fixed << std::setprecision(6) << hybrid_time
                      << std::setw(15) << std::fixed << std::setprecision(6) << direct_time
                      << std::setw(10) << std::fixed << std::setprecision(4) << omp_speedup
                      << std::setw(10) << std::fixed << std::setprecision(4) << mpi_speedup
                      << std::setw(10) << std::fixed << std::setprecision(4) << hybrid_speedup
                      << std::setw(10) << std::fixed << std::setprecision(4) << direct_speedup
                      << std::setw(15) << std::scientific << std::setprecision(3) << serial_l2 << "\n";

            if (n == 64)
                solver.save_vtk("solution_" + std::to_string(size) + "_n_" + std::to_string(n));
        }
    }

    // Only write results file on rank 0
    if (rank == 0)
    {
        std::ofstream ofs("test/data/results_" + std::to_string(size) + ".csv");
        ofs << "n,serial,omp,mpi,hybrid,direct,omp_speedup,mpi_speedup,hybrid_speedup,direct_speedup,l2_error\n";
        for (size_t i = 0; i < ns.size(); ++i)
        {
            ofs << ns[i] << "," << serial_times[i] << "," << omp_times[i] << "," << mpi_times[i] << "," << hybrid_times[i] << "," << direct_times[i] << ","
                << omp_speedups[i] << "," << mpi_speedups[i] << "," << hybrid_speedups[i] << "," << direct_speedups[i] << "," << l2_errors[i] << "\n";
        }
        ofs.close();
    }

    if (size == 4 && rank == 0)
    {
        std::cout << "========================" << std::endl;
        std::cout << "=== Plotting results ===" << std::endl;
        std::cout << "========================" << std::endl;
        plot::plot();
    }

    MPI_Finalize();
    return 0;
}