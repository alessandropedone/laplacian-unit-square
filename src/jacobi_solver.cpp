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

        // Update the boundary conditions

        // Check for convergence
        double residual = compute_error_serial(uh, previous);
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
    bool max_iter_reached = false;
    bool converged = false;

#ifdef _OPENMP
#pragma omp parallel num_threads(2) shared(uh, previous, max_iter_reached, converged)
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
                double residual = compute_error_serial(uh, previous);
                /* if (iteration % static_cast<int>(0.1 * max_iter) == 0)
                {
                    std::cout << std::setw(6) << iteration << std::setw(6) << "|| " << residual << std::endl;
                } */
                if (residual < tol)
                {
                    converged = true;
                    iter = ++iteration;
                }

                if (iteration == max_iter - 1)
                {
                    max_iter_reached = true;
                }
            }
        }
    }

    // Check if the maximum number of iterations was reached
    if (max_iter_reached)
    {
        std::cout << "Warning: Maximum number of iterations reached without convergence." << std::endl;
    }

    return;
}

void JacobiSolver::solve_mpi()
{
    // MPI implementation of the Jacobi solver
    // This function is a placeholder and should be implemented as needed
    //std::cout << "MPI implementation is not yet available." << std::endl;
    int initialized;
    MPI_Initialized(&initialized);

    if (initialized) {
        MPI_Comm mpi_comm = MPI_COMM_WORLD;

        int mpi_rank;
        MPI_Comm_rank(mpi_comm, &mpi_rank);

        int mpi_size;
        MPI_Comm_size(mpi_comm, &mpi_size);

        // To divide the work among processes
        unsigned int count = n / mpi_size;
        int remainder = n - count * mpi_size;

        // Vectors to store the number of elements to send to each
        // processor and the offset index where to start reading them from.
        std::vector<int> send_counts(mpi_size, 0);
        std::vector<int> send_start_idx(mpi_size, 0);

        // Vectors to store the number of elements to receive from each
        // processor and the offset index where to start writing them into.
        std::vector<int> recv_counts(mpi_size, 0);
        std::vector<int> recv_start_idx(mpi_size, 0);
        unsigned start_idx{0};

        // The first and last processors will have only one extra row,
        // while the others will have two ghost rows.
        if (mpi_rank == 0){

            recv_counts[0] = (0 < remainder) ? (count + 2) : count + 1;
            send_counts[0] = recv_counts[0] * n;

            recv_start_idx[0] = start_idx;
            send_start_idx[0] = start_idx * n;

            start_idx += recv_counts[0] - 1; //consider repetition of ghost row

            for (int i = 1; i < mpi_size - 1; ++i){

                recv_counts[i] = (i < remainder) ? (count + 3) : count + 2;
                send_counts[i] = recv_counts[i] * n;

                recv_start_idx[i] = start_idx;
                send_start_idx[i] = start_idx * n;
 
                start_idx += recv_counts[i] - 1; //consider repetition of ghost row
            }

            recv_counts[mpi_size-1] = (mpi_size-1 < remainder) ? (count + 2) : count + 1;
            send_counts[mpi_size-1] = recv_counts[mpi_size-1] * n;

            recv_start_idx[mpi_size-1] = start_idx;
            send_start_idx[mpi_size-1] = start_idx * n;

            start_idx += recv_counts[mpi_size-1] - 1; //consider repetition of ghost row

            for (int i = 0; i < mpi_size; ++i)
            {
                std::cout << "Rank " << i << " send_counts[" << i << "] = " << send_counts[i] << std::endl;
                std::cout << "Rank " << i << " recv_counts[" << i << "] = " << recv_counts[i] << std::endl;
                std::cout << "Rank " << i << " send_start_idx[" << i << "] = " << send_start_idx[i] << std::endl;
                std::cout << "Rank " << i << " recv_start_idx[" << i << "] = " << recv_start_idx[i] << std::endl;
            }
        }

        unsigned int local_rows = (mpi_rank < remainder) ? (count + 1) : count;
        local_rows += (mpi_rank == 0 || mpi_rank == mpi_size - 1) ? 1 : 2; // Add ghost rows
        //std::cout << "Number of rows on rank " << mpi_rank << ": " << local_rows << std::endl;
        MPI_Barrier(mpi_comm);

        /* int MPI_Scatterv(
                void *sendbuf,
                int *sendcnts,
                int *displs,
                MPI_Datatype sendtype,
                void *recvbuf,
                int recvcnt,
                MPI_Datatype recvtype,
                int root,
                MPI_Comm comm
        ); */
        std::vector<double> local_sol(local_rows * n);
        MPI_Scatterv(uh.data(),
               send_counts.data(),
               send_start_idx.data(),
               MPI_DOUBLE,
               local_sol.data(),
               local_rows * n,
               MPI_DOUBLE,
               0,
               mpi_comm);

        














        MPI_Gatherv(local_sol.data(),
                local_rows,
                MPI_DOUBLE,
                uh.data(),
                recv_counts.data(),
                recv_start_idx.data(),
                MPI_DOUBLE,
                0,
                mpi_comm);

        //if (mpi_rank == 0 && iter % 10 == 0)
          //  std::cout << "Iter " << iter << " Residual: " << globalResidual << std::endl;


        return;
    }
    else {
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

double JacobiSolver::compute_error_serial(const std::vector<double> &sol1, const std::vector<double> &sol2) const
{
    double error{0.0};
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            error += (sol1[i * n + j] - sol2[i * n + j]) * (sol1[i * n + j] - sol2[i * n + j]);
        }
    }
    error = std::sqrt(1.0 / (n - 1) * error);
    return error;
}

double JacobiSolver::compute_error_omp(const std::vector<double> &sol1, const std::vector<double> &sol2) const
{
    double error{0.0};
#ifdef _OPENMP
    int num_threads = omp_get_num_threads();
    int chunk_size = (n * n) / (num_threads); // ensures all elements are covered
#pragma omp parallel for collapse(2) schedule(static, chunk_size) reduction(+ : error) num_threads(num_threads)
#endif
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            error += (sol1[i * n + j] - sol2[i * n + j]) * (sol1[i * n + j] - sol2[i * n + j]);
        }
    }
    error = std::sqrt(1.0 / (n - 1) * error);
    return error;
}

double JacobiSolver::compute_error_serial(const std::vector<double> &sol1, const std::function<double(double, double)> &sol2) const
{
    double error{0.0};
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            error += (sol1[i * n + j] - fun_at(sol2, i, j)) * (sol1[i * n + j] - fun_at(sol2, i, j));
        }
    }
    error = std::sqrt(1.0 / (n - 1) * error);
    return error;
}

double JacobiSolver::compute_error_omp(const std::vector<double> &sol1, const std::function<double(double, double)> &sol2) const
{
    double error{0.0};
#ifdef _OPENMP
    int num_threads = omp_get_num_threads();
    int chunk_size = (n * n) / (num_threads); // ensures all elements are covered
#pragma omp parallel for collapse(2) schedule(static, chunk_size) reduction(+ : error) num_threads(num_threads)
#endif
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            error += (sol1[i * n + j] - fun_at(sol2, i, j)) * (sol1[i * n + j] - fun_at(sol2, i, j));
        }
    }
    error = std::sqrt(1.0 / (n - 1) * error);
    return error;
}
