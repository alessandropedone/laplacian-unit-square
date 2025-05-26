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
#include <mpi.h>

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
using Eigen::SparseMatrix;
using Eigen::Triplet;
using Eigen::VectorXd;


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

    // Initialize the converged variable
    bool converged = false;

    for (size_t iteration = 0; iteration < max_iter && !converged; ++iteration)
    {
        // Save the previous solution in a temporary vector
        std::copy(uh.begin(), uh.end(), previous.begin());

        // Perform the iteration
        for (size_t i = 1; i < n - 1; ++i)
        {
            for (size_t j = 1; j < n - 1; ++j)
            {
                // Jacobi iteration
                uh[i * n + j] = 0.25 * (previous[(i - 1) * n + j] + previous[(i + 1) * n + j] +
                                        previous[i * n + (j - 1)] + previous[i * n + (j + 1)] +
                                        h * h * fun_at(f, i, j));
            }
        }

        // Check for convergence
        double residual = compute_error_serial(uh, previous, n, n);
        if (residual < tol)
        {
            converged = true;
            // Save the number of iterations for convergence in iter memeber variable
            iter = ++iteration;
        }
        else if (iteration == max_iter - 1)
        {
            iter = ++iteration;
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

    // Initialize converged variable
    bool converged = false;

#ifdef _OPENMP
#pragma omp parallel num_threads(2) shared(uh, previous, converged)
#endif
    {

        for (size_t iteration = 0; iteration < max_iter && !converged; ++iteration)
        {
#ifdef _OPENMP
#pragma omp single
#endif
            {
                // Save the previous solution in a temporary vector
                std::copy(uh.begin(), uh.end(), previous.begin());
            }

#ifdef _OPENMP
            int num_threads = omp_get_num_threads();
            int chunk_size = (n * n) / (num_threads); // ensures all elements are covered
#pragma omp barrier
#pragma omp for collapse(2) schedule(static, chunk_size)
#endif
            // Perform the iteration
            for (size_t i = 1; i < n - 1; ++i)
            {
                for (size_t j = 1; j < n - 1; ++j)
                {
                    // Jacobi iteration
                    uh[i * n + j] = 0.25 * (previous[(i - 1) * n + j] + previous[(i + 1) * n + j] +
                                            previous[i * n + (j - 1)] + previous[i * n + (j + 1)] +
                                            h * h * fun_at(f, i, j));
                }
            }
#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
            {
                // Check for convergence
                double residual = compute_error_serial(uh, previous, n, n);
                if (residual < tol)
                {
                    converged = true;
                    // Save the number of iterations for convergence in iter memeber variable
                    iter = ++iteration;
                }
                else if (iteration == max_iter - 1)
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
    int initialized;
    MPI_Initialized(&initialized);

    if (initialized)
    {
        // We use MPI_COMM_WORLD ad communicator
        MPI_Comm mpi_comm = MPI_COMM_WORLD;

        // Get size and rank
        int mpi_rank, mpi_size;
        MPI_Comm_rank(mpi_comm, &mpi_rank);
        MPI_Comm_size(mpi_comm, &mpi_size);

        // Set the boundary conditions
        if(mpi_rank == 0){
            for (size_t i = 0; i < n; ++i)
            {
                uh[i] = fun_at(top_bc, 0, i);                       // Top boundary
                uh[i * n + (n - 1)] = fun_at(right_bc, i, n - 1);   // Right boundary
                uh[(n - 1) * n + i] = fun_at(bottom_bc, n - 1, i);  // Bottom boundary
                uh[i * n] = fun_at(left_bc, i, 0);                  // Left boundary
            }
        }

        // Compute these two quantities to divide the work among processes
        unsigned int count = n / mpi_size;
        int remainder = n - count * mpi_size;

        // Vector to store the number of elements to send to each processor,
        // and to store the number of elements to receive from each processor.
        std::vector<int> counts(mpi_size, 0);

        // Vector to store the offset index where to start reading/writing them from.
        std::vector<int> start_idxs(mpi_size, 0);

        // Auxiliary index
        unsigned start_idx{0};

        if (mpi_rank == 0)
        {
            if (mpi_size > 1)
            {
                // The first and last processors will have only one extra row,
                // while the others will have two ghost rows.

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
                // We handle the case with only one process separately
                counts[0] = n * n; // Only one process, all data
                start_idxs[0] = 0;
            }
        }

        // Compute the number of rows of the local grid
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

        // Broadcast start_idxs because the y will be used during computation
        MPI_Bcast(start_idxs.data(), mpi_size, MPI_INT, 0, mpi_comm);

        // Synchronize all processes
        MPI_Barrier(mpi_comm);

        // Declare the local grid
        std::vector<double> local_uh(local_rows * n);

        // Scatter the initial guess between processes
        MPI_Scatterv(uh.data(),
                     counts.data(),
                     start_idxs.data(),
                     MPI_DOUBLE,
                     local_uh.data(),
                     local_rows * n,
                     MPI_DOUBLE,
                     0,
                     mpi_comm);

        // Grid that will containt the solution at the previous iteration
        std::vector<double> local_previous(local_rows * n);

        // Define converged variable
        bool converged = false;

        // Define h
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
                    // Jacobi iteration
                    local_uh[i * n + j] = 0.25 * (local_previous[(i - 1) * n + j] + local_previous[(i + 1) * n + j] +
                                                local_previous[i * n + (j - 1)] + local_previous[i * n + (j + 1)] +
                                                h * h * fun_at(f, start_idxs[mpi_rank] / n + i, j));
                }
            }

            // Check for convergence
            // Compute the local residual
            double local_residual = compute_error_serial(local_uh, local_previous, local_rows, n);
            double global_residual;
            // Ensure all processes have computed their local residual before reduction
            MPI_Barrier(mpi_comm);
            // Find the maximum residual across all processes
            MPI_Allreduce(&local_residual, &global_residual, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
            // The method converged if all local residual satisfy the convergence criterion
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

            // Bidirectional ghost cell exchange
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

        // Gather the results from local grids in uh (global grid)
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

#ifndef _OPENMP
    std::cout << "Warning: OpenMP is not enabled. Falling back to serial execution." << std::endl;
#endif
    int initialized;
    MPI_Initialized(&initialized);

    if (initialized)
    {
        // We use MPI_COMM_WORLD ad communicator
        MPI_Comm mpi_comm = MPI_COMM_WORLD;

        // Get size and rank
        int mpi_rank, mpi_size;
        MPI_Comm_rank(mpi_comm, &mpi_rank);
        MPI_Comm_size(mpi_comm, &mpi_size);

        // Compute these two quantities to divide the work among processes
        unsigned int count = n / mpi_size;
        int remainder = n - count * mpi_size;

        // Vector to store the number of elements to send to each processor,
        // and to store the number of elements to receive from each processor.
        std::vector<int> counts(mpi_size, 0);

        // Vector to store the offset index where to start reading/writing them from.
        std::vector<int> start_idxs(mpi_size, 0);

        // Auxiliary index
        unsigned start_idx{0};

        if (mpi_rank == 0)
        {
            if (mpi_size > 1)
            {
                // The first and last processors will have only one extra row,
                // while the others will have two ghost rows.

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
                // We handle the case with only one process separately
                counts[0] = n * n; // Only one process, all data
                start_idxs[0] = 0;
            }
        }

        // Compute the number of rows of the local grid
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

        // Broadcast start_idxs because the y will be used during computation
        MPI_Bcast(start_idxs.data(), mpi_size, MPI_INT, 0, mpi_comm);

        // Synchronize all processes
        MPI_Barrier(mpi_comm);

        // Declare the local grid
        std::vector<double> local_uh(local_rows * n);

        // Scatter the initial guess between processes
        MPI_Scatterv(uh.data(),
                     counts.data(),
                     start_idxs.data(),
                     MPI_DOUBLE,
                     local_uh.data(),
                     local_rows * n,
                     MPI_DOUBLE,
                     0,
                     mpi_comm);

        // Grid that will containt the solution at the previous iteration
        std::vector<double> local_previous(local_rows * n);

        // Define converged variable
        bool converged = false;

        // Define h
        const double h = 1.0 / (n - 1);

#ifdef _OPENMP
#pragma omp parallel num_threads(2) shared(local_uh, local_previous, converged)
#endif
        for (size_t iteration = 0; iteration < max_iter && !converged; ++iteration)
        {
#ifdef _OPENMP
#pragma omp single
#endif
            {
                // Save the previous solution for convergence check
                std::copy(local_uh.begin(), local_uh.end(), local_previous.begin());
            }
#ifdef _OPENMP
            int num_threads = omp_get_num_threads();
            int chunk_size = (local_rows * n) / (num_threads); // ensures all elements are covered
#pragma omp barrier
#pragma omp for collapse(2) schedule(static, chunk_size)
#endif
            // Perform the iteration
            for (size_t i = 1; i < local_rows - 1; ++i)
            {
                for (size_t j = 1; j < n - 1; ++j)
                {
                    // Jacobi iteration
                    local_uh[i * n + j] = 0.25 * (local_previous[(i - 1) * n + j] + local_previous[(i + 1) * n + j] +
                                                  local_previous[i * n + (j - 1)] + local_previous[i * n + (j + 1)] +
                                                  h * h * fun_at(f, start_idxs[mpi_rank] / n + i, j));
                }
            }
#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
            {
                // Check for convergence
                // Compute the local residual
                double local_residual = compute_error_serial(local_uh, local_previous, local_rows, n);
                double global_residual;
                // Ensure all processes have computed their local residual before reduction
                MPI_Barrier(mpi_comm);
                // Find the maximum residual across all processes
                MPI_Allreduce(&local_residual, &global_residual, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
                // The method converged if all local residual satisfy the convergence criterion
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

                // Bidirectional ghost cell exchange
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
        }

        // Synchronize all processes before gathering results
        MPI_Barrier(mpi_comm);

        // Gather the results from local grids in uh (global grid)
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

void JacobiSolver::solve_direct()
{
    int initialized;
    MPI_Initialized(&initialized);

    if (initialized)
    {
        // We use MPI_COMM_WORLD ad communicator
        MPI_Comm mpi_comm = MPI_COMM_WORLD;

        // Get size and rank
        int mpi_rank, mpi_size;
        MPI_Comm_rank(mpi_comm, &mpi_rank);
        MPI_Comm_size(mpi_comm, &mpi_size);

        // Set the boundary conditions
        if(mpi_rank == 0){
            for (size_t i = 0; i < n; ++i)
            {
                uh[i] = fun_at(top_bc, 0, i);                       // Top boundary
                uh[i * n + (n - 1)] = fun_at(right_bc, i, n - 1);   // Right boundary
                uh[(n - 1) * n + i] = fun_at(bottom_bc, n - 1, i);  // Bottom boundary
                uh[i * n] = fun_at(left_bc, i, 0);                  // Left boundary
            }
        }

        // Compute these two quantities to divide the work among processes
        unsigned int count = n / mpi_size;
        int remainder = n - count * mpi_size;

        // Vector to store the number of elements to send to each processor,
        // and to store the number of elements to receive from each processor.
        std::vector<int> counts(mpi_size, 0);

        // Vector to store the offset index where to start reading/writing them from.
        std::vector<int> start_idxs(mpi_size, 0);

        // Auxiliary index
        unsigned start_idx{0};

        if (mpi_rank == 0)
        {
            if (mpi_size > 1)
            {
                // The first and last processors will have only one extra row,
                // while the others will have two ghost rows.

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
                // We handle the case with only one process separately
                counts[0] = n * n; // Only one process, all data
                start_idxs[0] = 0;
            }
        }

        // Compute the number of rows of the local grid
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

        // Broadcast start_idxs because the y will be used during computation
        MPI_Bcast(start_idxs.data(), mpi_size, MPI_INT, 0, mpi_comm);

        // Synchronize all processes
        MPI_Barrier(mpi_comm);

        // Declare the local grid
        std::vector<double> local_uh(local_rows * n);

        // Scatter the initial guess between processes
        MPI_Scatterv(uh.data(),
                     counts.data(),
                     start_idxs.data(),
                     MPI_DOUBLE,
                     local_uh.data(),
                     local_rows * n,
                     MPI_DOUBLE,
                     0,
                     mpi_comm);

        // Grid that will containt the solution at the previous iteration
        std::vector<double> local_previous(local_rows * n);

        // Define converged variable
        bool converged = false;

        // Define h
        const double h = 1.0 / (n - 1);
   
        for (size_t iteration = 0; iteration < max_iter && !converged; ++iteration)
        {
            // Save the previous solution for convergence check
            std::copy(local_uh.begin(), local_uh.end(), local_previous.begin());

            // Assemble local system
            unsigned working_rows = local_rows - 2; // in each case, top and bottom rows are given
            unsigned working_cols = n - 2; // in each case, left and right columns are given
            // Unknowns are in the interior of the grid, so we have working_rows * working_cols unknowns
            SparseMatrix<double> A(working_rows * working_cols, working_rows * working_cols);
            std::vector<Triplet<double>> triplets; // to insert elements in sparse matrix
            VectorXd b = VectorXd::Zero(working_rows * working_cols);

            for (size_t i = 0; i < working_rows; ++i)
            {
                for (size_t j = 0; j < working_cols; ++j)
                {
                    int idx = i * working_cols + j;
                    triplets.emplace_back(idx, idx, 4.0);
                    if (i > 0) {// Use top neighbor
                        triplets.emplace_back(idx, idx - working_cols, -1.0); 
                    }
                    else {
                        // If we are at the first local row, use upper ghost row
                        b(idx) += local_uh[j]; 
                    }
                    if (i < working_rows - 1) { // Use bottom neighbor
                        triplets.emplace_back(idx, idx + working_cols, -1.0);
                    }
                    else {
                        // If we are at the last local row, use lower ghost row
                        b(idx) += local_uh[(local_rows - 1) * n + j];
                    }
                    if (j > 0) { // Use left neighbor
                        triplets.emplace_back(idx, idx - 1, -1.0);
                    }
                    else {
                        // If we are at the first local column, use left ghost column
                        b(idx) += local_uh[i * n];
                    }
                    if (j < working_cols - 1) { // Use right neighbor
                        triplets.emplace_back(idx, idx + 1, -1.0);
                    }
                    else {
                        // If we are at the last local column, use right ghost column
                        b(idx) += local_uh[i * n + (n - 1)];
                    }
                    b(idx) += h * h * fun_at(f, start_idxs[mpi_rank] / n + i + 1, j + 1);
                }
            }

            A.setFromTriplets(triplets.begin(), triplets.end());

            // Solve the local system
            Eigen::SimplicialLLT<SparseMatrix<double>> solver;
            solver.compute(A);
            VectorXd x = solver.solve(b);

            for (unsigned i = 1; i < local_rows - 1; ++i)
                for (unsigned j = 1; j < n - 1; ++j)
                    local_uh[i * n + j] = x[(i - 1) * working_cols + j];
        
            // Check for convergence
            // Compute the local residual
            double local_residual = compute_error_serial(local_uh, local_previous, local_rows, n);
            double global_residual;
            // Ensure all processes have computed their local residual before reduction
            MPI_Barrier(mpi_comm);
            // Find the maximum residual across all processes
            MPI_Allreduce(&local_residual, &global_residual, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
            // The method converged if all local residual satisfy the convergence criterion
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

            // Bidirectional ghost cell exchange
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

        // Gather the results from local grids in uh (global grid)
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

double JacobiSolver::compute_error_serial(const std::vector<double> &sol1, const std::function<double(std::vector<double>)> &sol2, unsigned rows, unsigned cols) const
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

double JacobiSolver::compute_error_omp(const std::vector<double> &sol1, const std::function<double(std::vector<double>)> &sol2, unsigned rows, unsigned cols) const
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
