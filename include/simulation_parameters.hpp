

/**
 * @file simulation_parameters.hpp
 * @brief Header file containing simulation parameters structure and string broadcasting utility for MPI-based Poisson equation solver.
 *
 * This file defines the SimulationParameters structure that holds all configuration values
 * needed for solving the Poisson equation, including boundary conditions, solver parameters,
 * and mathematical functions. It also provides MPI broadcasting capabilities to distribute
 * these parameters across all processes in a parallel computation environment.
 *
 * The parameters include:
 * - Right-hand side function and exact solution (as muParserX expressions)
 * - Boundary conditions for all four edges of the computational domain
 * - Numerical solver settings (tolerance and maximum iterations)
 */
#ifndef SIMULATION_PARAMETERS_HPP
#define SIMULATION_PARAMETERS_HPP

#include <string>
#include <mpi.h>

/**
 * @namespace solver
 * @brief Contains utilities and data structures for distributed numerical simulations.
 * 
 * This namespace provides MPI-based broadcasting functionality and parameter management
 * for parallel Poisson equation solvers. It includes utilities for distributing
 * simulation configuration across multiple processes in a distributed computing environment.
 */
namespace solver
{
    /// @brief  Broadcast a string across all MPI processes.
    /// @param str The string to broadcast.
    /// @param root The rank of the root process that provides the string.
    /// @param comm The MPI communicator to use for broadcasting.
    void broadcast_string(std::string &str, int root, MPI_Comm comm)
    {
        int rank;
        MPI_Comm_rank(comm, &rank);

        int str_len = (rank == root) ? str.length() : 0;
        MPI_Bcast(&str_len, 1, MPI_INT, root, comm);

        if (rank != root)
        {
            str.resize(str_len);
        }

        if (str_len > 0)
        {
            MPI_Bcast(&str[0], str_len, MPI_CHAR, root, comm);
        }
    }

    /// @brief  Simulation parameters structure for holding configuration values.
    /// @details This structure holds the parameters read from the data file, including
    ///          the right-hand side function, exact solution, boundary conditions,
    ///          tolerance, and maximum iterations. It provides a method to broadcast
    ///          these parameters to all MPI processes.
    struct SimulationParameters
    {
        /// @brief Right-hand side function as a string.
        /// @details This string represents the function f(x, y) used in the Poisson equation.
        ///          It is expected to be a valid muParserX expression.
        std::string f_str;

        /// @brief Exact solution function as a string.
        /// @details This string represents the exact solution u(x, y) of the Poisson equation.
        std::string uex_str;

        /// @brief Boundary condition for the top edge as a string.
        std::string bc_top_str;

        /// @brief Boundary condition for the right edge as a string.
        std::string bc_right_str;

        /// @brief Boundary condition for the bottom edge as a string.
        std::string bc_bottom_str;

        /// @brief Boundary condition for the left edge as a string.
        std::string bc_left_str;

        /// @brief Tolerance for convergence.
        double tol;

        /// @brief Maximum number of iterations for the Jacobi solver.
        unsigned max_iter;

        /// @brief Broadcast the simulation parameters to all MPI processes.
        /// @param root The rank of the root process that provides the parameters.
        /// @param comm  The MPI communicator to use for broadcasting.
        void broadcast(int root, MPI_Comm comm)
        {
            broadcast_string(f_str, root, comm);
            broadcast_string(uex_str, root, comm);
            broadcast_string(bc_top_str, root, comm);
            broadcast_string(bc_right_str, root, comm);
            broadcast_string(bc_bottom_str, root, comm);
            broadcast_string(bc_left_str, root, comm);

            MPI_Bcast(&tol, 1, MPI_DOUBLE, root, comm);
            MPI_Bcast(&max_iter, 1, MPI_UNSIGNED, root, comm);
        }
    };
}
#endif // SIMULATION_PARAMETERS_HPP