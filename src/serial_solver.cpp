/// @file serial_solver.cpp
/// @brief This file contains the implementation of the SerialSolver class.
/// @details The class implements an iterative solver for a given equation.
/// @details The class provides methods to set the boundary conditions,
///          initial guess, exact solution, and right-hand side of the equation.    
///          It also includes a method to print the computed solution.


#include <iostream>
#include <vector>
#include <cmath>
#include "serial_solver.hpp"

SerialSolver::SerialSolver(const std::vector<double>& exact_sol,
                    const std::vector<double>& initial_guess,
                    const std::vector<double>& rhs,
                    const std::vector<double>& topbc,
                    const std::vector<double>& rightbc,
                    const std::vector<double>& bottombc,
                    const std::vector<double>& leftbc,
                    size_t n, unsigned max_iter, double tol):
                            exact_sol(exact_sol),
                            sol(initial_guess),
                            rhs(rhs),
                            topbc(topbc),
                            rightbc(rightbc),
                            bottombc(bottombc),
                            leftbc(leftbc),
                            n(n),
                            max_iter(max_iter),
                            tol(tol) {};

SerialSolver::~SerialSolver() {};

void SerialSolver::solve(const std::vector<double>& x_points, const std::vector<double>& y_points) {
    std::cout << "Solving the equation iteratively..." << std::endl;
    // Initialize the grid size
    const double h = 1.0 / (n - 1);

    // Set the boundary conditions
    for (size_t i = 0; i < n; ++i) {
            sol[i] = topbc[i]; // Top boundary
            sol[(n - 1) * n + i] = bottombc[i]; // Bottom boundary
            sol[i * n] = leftbc[i]; // Left boundary
            sol[i * n + (n - 1)] = rightbc[i]; // Right boundary
        }
    for (n_iter = 0; n_iter < max_iter; ++n_iter) {
        // Save the previous solution for convergence check
        std::vector<double> previous(sol);
        // Perform the iteration
        for (size_t i = 1; i < n - 1; ++i) {
            for (size_t j = 1; j < n - 1; ++j) {
                sol[i * n + j] = 0.25 * (previous[(i - 1) * n + j] + previous[(i + 1) * n + j] +
                                         previous[i * n + (j - 1)] + previous[i * n + (j + 1)] +
                                         h * h * rhs[i * n + j]);
            }
        }
        // Check for convergence
        double residual = compute_error(h, previous);
        if (residual < tol) {
            break;
        }
    }
    if (n_iter == max_iter) {
        std::cout << "Warning: Maximum number of iterations reached without convergence." << std::endl;
    } else {
        std::cout << "Converged in " << n_iter << " iterations." << std::endl;
    }
    // Print the computed solution
    std::cout << "Computed solution:" << std::endl;
    SerialSolver::printSolution();
    // Compute the error
    double error = compute_error(h, exact_sol);
    std::cout << "Error: " << error << std::endl;
    return;
};

void SerialSolver::set_bc(const std::vector<double>& topbc,
                          const std::vector<double>& rightbc,
                          const std::vector<double>& bottombc,
                          const std::vector<double>& leftbc) {

    this->topbc = topbc;
    this->rightbc = rightbc;
    this->bottombc = bottombc;
    this->leftbc = leftbc;
};

void SerialSolver::set_initial_guess(const std::vector<double>& initial_guess){
    this->sol = initial_guess;
};

void SerialSolver::set_exact_sol(const std::vector<double>& exact_sol){
    this->exact_sol = exact_sol;
};

void SerialSolver::set_rhs(const std::vector<double>& rhs){
    this->rhs = rhs;
};

double SerialSolver::compute_error(const double h, const std::vector<double> & reference) const{
    double error{0.0};
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            error += (sol[i * n + j] - reference[i * n + j])*(sol[i * n + j] - reference[i * n + j]);
        }
    }
    error = std::sqrt(h * error);   
    return error;
};

void SerialSolver::printSolution() const{
    std::cout << "Computed solution:" << std::endl;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            std::cout << sol[i * n + j] << "\t ";
        }
        std::cout << std::endl;
    }
};
