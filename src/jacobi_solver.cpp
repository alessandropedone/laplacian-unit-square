/// @file jacobi_serial_solver.cpp
/// @brief This file contains the implementation of the JacobiSolver class.
/// @details The class implements an iterative solver for a given equation.
/// @details The class provides methods to set the boundary conditions,
///          initial guess, exact solution, and right-hand side of the equation.
///          It also includes a method to print the computed solution.

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>

#include "jacobi_solver.hpp"

void JacobiSolver::solve_serial()
{

    // Initialize h
    const double h = 1.0 / (n - 1);

    // Set the boundary conditions
    for (size_t i = 0; i < n; ++i)
    {
        uh[i] = fun_at(top_bc, i, n - 1);                 // Top boundary
        uh[i * n + (n - 1)] = fun_at(right_bc, n - 1, i); // Right boundary
        uh[(n - 1) * n + i] = fun_at(bottom_bc, i, 0);    // Bottom boundary
        uh[i * n] = fun_at(left_bc, 0, i);                // Left boundary
    }

    // Initialize the previous solution vector
    std::vector<double> previous(n * n);
    // std::cout << std::setw(6) << "Iteration" << std::setw(6) << "|| Residual" << std::endl;
    bool converged = false;

    for (size_t iteration = 0; iteration < max_iter && !converged; ++iteration)
    {
        // Save the previous solution for convergence check
        std::copy(uh.begin(), uh.end(), previous.begin());

        // Perform the iteration
        for (size_t i = 1; i < n - 1; ++i)
        {
            for (size_t j = 1; j < n - 1; ++j)
            {
                uh[i * n + j] = 0.25 * (previous[(i - 1) * n + j] + previous[(i + 1) * n + j] +
                                        previous[i * n + (j - 1)] + previous[i * n + (j + 1)] +
                                        h * h * fun_at(f, i, j));
            }
        }

        // Check for convergence
        double residual = compute_error_serial(uh, previous, n, n);
        /* if (iteration % static_cast<int>(0.1 * max_iter) == 0)
        {
            std::cout << std::setw(6) << iteration << std::setw(6) << "|| " << residual << std::endl;
        }  */
        if (residual < tol)
        {
            converged = true;
            iter = ++iteration;
        }
        if (iteration == max_iter - 1)
        {
            std::cout << "Warning: Maximum number of iterations reached without convergence." << std::endl;
        }
    }
    return;
}

void JacobiSolver::solve_omp()
{

#ifndef _OPENMP
    std::cout << "Warning: OpenMP is not enabled. Falling back to serial execution." << std::endl;
#endif

    // Initialize h
    const double h = 1.0 / (n - 1);

    // Set the boundary conditions
    for (size_t i = 0; i < n; ++i)
    {
        uh[i] = fun_at(top_bc, i, n - 1);                 // Top boundary
        uh[i * n + (n - 1)] = fun_at(right_bc, n - 1, i); // Right boundary
        uh[(n - 1) * n + i] = fun_at(bottom_bc, i, 0);    // Bottom boundary
        uh[i * n] = fun_at(left_bc, 0, i);                // Left boundary
    }

    // Initialize the previous solution vector
    std::vector<double> previous(n * n);
    // std::cout << std::setw(6) << "Iteration" << std::setw(6) << "|| Residual" << std::endl;
    bool converged = false;

#ifdef _OPENMP
#pragma omp parallel num_threads(2) shared(uh, previous, converged)
#endif
    {
        // Initialize the previous solution vector
        for (size_t iteration = 0; iteration < max_iter && !converged; ++iteration)
        {
// Save the previous solution for convergence check
#ifdef _OPENMP
#pragma omp single
#endif
            {
                std::copy(uh.begin(), uh.end(), previous.begin());
            }
            // Perform the iteration
#ifdef _OPENMP
            int num_threads = omp_get_num_threads();
            int chunk_size = (n * n) / (num_threads); // ensures all elements are covered
#pragma omp barrier

            // Update the solution using the Jacobi method
#pragma omp for collapse(2) schedule(static, chunk_size)
#endif
            for (size_t i = 1; i < n - 1; ++i)
            {
                for (size_t j = 1; j < n - 1; ++j)
                {
                    uh[i * n + j] = 0.25 * (previous[(i - 1) * n + j] + previous[(i + 1) * n + j] +
                                            previous[i * n + (j - 1)] + previous[i * n + (j + 1)] +
                                            h * h * fun_at(f, i, j));
                }
            }
#ifdef _OPENMP
#pragma omp barrier

            // Update the boundary conditions
// Check for convergence
#pragma omp single
#endif
            {
                double residual = compute_error_serial(uh, previous, n, n);
                /* if (iteration % static_cast<int>(0.1 * max_iter) == 0)
                {
                    std::cout << std::setw(6) << iteration << std::setw(6) << "|| " << residual << std::endl;
                } */
                if (residual < tol)
                {
                    converged = true;
                    iter = ++iteration;
                }
                else if (iteration == max_iter - 1) // Check if the maximum number of iterations was reached
                {
                    iter = ++iteration;
                    std::cout << "Warning: Maximum number of iterations reached without convergence." << std::endl;
                }
            }
        }
    }

    return;
}

void JacobiSolver::solve_mpi()
{
    // MPI implementation of the Jacobi solver
    // This function is a placeholder and should be implemented as needed
    // std::cout << "MPI implementation is not yet available." << std::endl;
    int initialized;
    MPI_Initialized(&initialized);

    if (initialized)
    {
        MPI_Comm mpi_comm = MPI_COMM_WORLD;

        int mpi_rank;
        MPI_Comm_rank(mpi_comm, &mpi_rank);

        int mpi_size;
        MPI_Comm_size(mpi_comm, &mpi_size);

        // To divide the work among processes
        unsigned int count = n / mpi_size;
        int remainder = n - count * mpi_size;

        // Vector to store the number of elements to send to each processor,
        // and to store the number of elements to receive from each processor.
        std::vector<int> counts(mpi_size, 0);
        // Vector to store the offset index where to start reading/writing them from.
        std::vector<int> start_idxs(mpi_size, 0);

        unsigned start_idx{0};

        // The first and last processors will have only one extra row,
        // while the others will have two ghost rows.
        if (mpi_rank == 0)
        {
            if (mpi_size > 1)
            {
                counts[0] = ((0 < remainder) ? (count + 1 + 1) : count + 1) * n;
                start_idxs[0] = 0;
                start_idx = counts[0] - 2 * n;

                for (int i = 1; i < mpi_size - 1; ++i)
                {
                    counts[i] = ((i < remainder) ? (count + 1 + 2) : count + 2) * n;
                    start_idxs[i] = start_idx;
                    start_idx += counts[i] - 2 * n; // consider repetition of ghost row
                }
                counts[mpi_size - 1] = ((mpi_size - 1 < remainder) ? (count + 1 + 1) : count + 1) * n;
                start_idxs[mpi_size - 1] = start_idx;
            }
            else
            {
                counts[0] = n * n; // Only one process, all data
                start_idxs[0] = 0;
            }

            for (int i = 0; i < mpi_size; ++i)
            {
                std::cout << "Rank " << i << " start_idxs[" << i << "] = " << start_idxs[i] << std::endl;
                std::cout << "Rank " << i << " counts[" << i << "] = " << counts[i] << std::endl;
            }
        }

        unsigned int local_rows;
        if (mpi_size > 1)
        {
            local_rows = (mpi_rank < remainder) ? (count + 1) : count;
            local_rows += (mpi_rank == 0 || mpi_rank == mpi_size - 1) ? 1 : 2; // Add ghost rows
        }
        else
        {
            local_rows = n; // Only one process, all rows
        }

        MPI_Bcast(start_idxs.data(), mpi_size, MPI_INT, 0, mpi_comm);
        MPI_Barrier(mpi_comm);

        std::vector<double> local_uh(local_rows * n);
        MPI_Scatterv(uh.data(),
                     counts.data(),
                     start_idxs.data(),
                     MPI_DOUBLE,
                     local_uh.data(),
                     local_rows * n,
                     MPI_DOUBLE,
                     0,
                     mpi_comm);

        std::vector<double> local_previous(local_rows * n);
        bool converged = false;
        const double h = 1.0 / (n - 1);

        for (size_t iteration = 0; iteration < max_iter && !converged; ++iteration)
        {
            // Save the previous solution for convergence check
            std::copy(local_uh.begin(), local_uh.end(), local_previous.begin());

            // Perform the iteration
            for (size_t i = 1; i < local_rows - 1; ++i)
            {
                for (size_t j = 1; j < n - 1; ++j)
                {
                    local_uh[i * n + j] = 0.25 * (local_previous[(i - 1) * n + j] + local_previous[(i + 1) * n + j] +
                                                  local_previous[i * n + (j - 1)] + local_previous[i * n + (j + 1)] +
                                                  h * h * fun_at(f, start_idxs[mpi_rank] / n + i, j));
                }
            }

            // Check for convergence
            double local_residual = compute_error_serial(local_uh, local_previous, local_rows, n);
            double global_residual;
            // Ensure all processes have computed their local residual before reduction
            MPI_Barrier(mpi_comm);
            // Find the maximum residual across all processes
            MPI_Allreduce(&local_residual, &global_residual, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
            converged = (global_residual < tol);
            if (converged)
            {
                iter = ++iteration;
            }
            else if (iteration == max_iter - 1)
            {
                iter = ++iteration;
                std::cout << "Warning: Maximum number of iterations reached without convergence." << std::endl;
            }

            // Proper bidirectional ghost cell exchange
            if (mpi_size > 1)
            {
                // Send/receive with next rank
                if (mpi_rank < mpi_size - 1)
                {
                    // Send my last interior row to next rank's ghost row
                    MPI_Send(&local_uh[(local_rows - 2) * n], n, MPI_DOUBLE, mpi_rank + 1, 0, mpi_comm);
                    // Receive next rank's first interior row into my ghost row
                    MPI_Recv(&local_uh[(local_rows - 1) * n], n, MPI_DOUBLE, mpi_rank + 1, 0, mpi_comm, MPI_STATUS_IGNORE);
                }

                // Send/receive with previous rank
                if (mpi_rank > 0)
                {
                    // Send my first interior row to previous rank's ghost row
                    MPI_Send(&local_uh[1 * n], n, MPI_DOUBLE, mpi_rank - 1, 0, mpi_comm);
                    // Receive previous rank's last interior row into my ghost row
                    MPI_Recv(&local_uh[0], n, MPI_DOUBLE, mpi_rank - 1, 0, mpi_comm, MPI_STATUS_IGNORE);
                }
            }
        }

        // Synchronize all processes before gathering results
        MPI_Barrier(mpi_comm);

        MPI_Gatherv(local_uh.data(),
                    local_rows * n,
                    MPI_DOUBLE,
                    uh.data(),
                    counts.data(),
                    start_idxs.data(),
                    MPI_DOUBLE,
                    0,
                    mpi_comm);
    }
    else
    {
        std::cerr << "Error: MPI is not initialized." << std::endl;
        return;
    }
}

void JacobiSolver::solve_hybrid()
{
    // Hybrid MPI and OpenMP implementation of the Jacobi solver
    // This function is a placeholder and should be implemented as needed
    std::cout << "Hybrid MPI and OpenMP implementation is not yet available." << std::endl;
}

double JacobiSolver::compute_error_serial(const std::vector<double> &sol1, const std::vector<double> &sol2, unsigned rows, unsigned cols) const
{
    double error{0.0};
    for (unsigned i = 0; i < rows; ++i)
    {
        for (unsigned j = 0; j < cols; ++j)
        {
            error += (sol1[i * n + j] - sol2[i * n + j]) * (sol1[i * n + j] - sol2[i * n + j]);
        }
    }
    error = std::sqrt(1.0 / (n - 1) * error);
    return error;
}

double JacobiSolver::compute_error_omp(const std::vector<double> &sol1, const std::vector<double> &sol2, unsigned rows, unsigned cols) const
{
    double error{0.0};
#ifdef _OPENMP
    int num_threads = omp_get_num_threads();
    int chunk_size = (n * n) / (num_threads); // ensures all elements are covered
#pragma omp parallel for collapse(2) schedule(static, chunk_size) reduction(+ : error) num_threads(num_threads)
#endif
    for (unsigned i = 0; i < rows; ++i)
    {
        for (unsigned j = 0; j < cols; ++j)
        {
            error += (sol1[i * n + j] - sol2[i * n + j]) * (sol1[i * n + j] - sol2[i * n + j]);
        }
    }
    error = std::sqrt(1.0 / (n - 1) * error);
    return error;
}

double JacobiSolver::compute_error_serial(const std::vector<double> &sol1, const std::function<double(double, double)> &sol2, unsigned rows, unsigned cols) const
{
    double error{0.0};
    for (unsigned i = 0; i < rows; ++i)
    {
        for (unsigned j = 0; j < cols; ++j)
        {
            error += (sol1[i * n + j] - fun_at(sol2, i, j)) * (sol1[i * n + j] - fun_at(sol2, i, j));
        }
    }
    error = std::sqrt(1.0 / (n - 1) * error);
    return error;
}

double JacobiSolver::compute_error_omp(const std::vector<double> &sol1, const std::function<double(double, double)> &sol2, unsigned rows, unsigned cols) const
{
    double error{0.0};
#ifdef _OPENMP
    int num_threads = omp_get_num_threads();
    int chunk_size = (n * n) / (num_threads); // ensures all elements are covered
#pragma omp parallel for collapse(2) schedule(static, chunk_size) reduction(+ : error) num_threads(num_threads)
#endif
    for (unsigned i = 0; i < rows; ++i)
    {
        for (unsigned j = 0; j < cols; ++j)
        {
            error += (sol1[i * n + j] - fun_at(sol2, i, j)) * (sol1[i * n + j] - fun_at(sol2, i, j));
        }
    }
    error = std::sqrt(1.0 / (n - 1) * error);
    return error;
}
