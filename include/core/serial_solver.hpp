/// @file serial_solver.hpp
/// @brief Header file for the SerialSolver class
/// @date 2025-18-05
/// @details This file contains the declaration of the SerialSolver class,
///          which implements an iterative solver for a given equation.
///          The class provides methods to set the boundary conditions,
///          initial guess, exact solution, and right-hand side of the equation.        
///          It also includes a method to print the computed solution.
/// @details The class is designed to be used in a serial computing environment.
///          It is not optimized for parallel computing and does not use any
///          parallel computing libraries or techniques.
/// @details The class is intended to be used as a comparison for more complex
///          solvers that may incorporate parallelism or other advanced features.
#ifndef     SERIAL_SOLVER_HPP
#define     SERIAL_SOLVER_HPP
#include <iostream>
#include <vector>

class SerialSolver {
    public:
        SerialSolver() = default;
        SerialSolver(const std::vector<double>& exact_sol,
                    const std::vector<double>& initial_guess,
                    const std::vector<double>& rhs,
                    const std::vector<double>& topbc,
                    const std::vector<double>& rightbc,
                    const std::vector<double>& bottombc,
                    const std::vector<double>& leftbc,
                    size_t n, unsigned max_iter, double tol);

        /// @brief default destructor
        ~SerialSolver();
        

        /// @brief implement iterative solver for the equation
        /// @param x_points points of discretization in x direction
        /// @param y_points points of discretization in y direction
        void solve(const std::vector<double>& x_points, const std::vector<double>& y_points);

        // SETTERS

        /// @brief set the boundary conditions of the system
        /// @param bc boundary conditions
        void set_bc(const std::vector<double>& topbc,
                    const std::vector<double>& rightbc,
                    const std::vector<double>& bottombc,
                    const std::vector<double>& leftbc);

        /// @brief set the initial guess for the solution
        /// @param initial_guess initial guess for the solution
        void set_initial_guess(const std::vector<double>& initial_guess);

        /// @brief set the exact solution of the equation
        /// @param exact_sol exact solution of the equation
        void set_exact_sol(const std::vector<double>& exact_sol);

        /// @brief set the right-hand side of the equation  
        /// @param rhs right-hand side of the equation
        void set_rhs(const std::vector<double>& rhs);
        

        /// @brief set the number of max iterations
        /// @param max_iter maximum number of iterations
        void set_max_iter(size_t max_iter) {
            this->max_iter = max_iter;
        };
        /// @brief set the tolerance for convergence
        /// @param tol tolerance for convergence
        void set_tol(double tol) {
            this->tol = tol;
        };

        // GETTERS

        /// @brief get the number of iterations
        /// @return number of iterations
        unsigned get_n_iter() const {
            return n_iter;
        }
        
        /// @brief get the computed solution
        /// @return computed solution
        const std::vector<double> get_solution() const {
            return sol;
        };

        /// @brief compute the error of the computed solution
        /// @param h grid size
        /// @param reference function to compare the computed solution with
        /// @return error of the computed solution
        /// @note It can be used to compare the computed solution with the exact solution
        ///       but also with the previous approximated solution in solve method
        double compute_error(const double h, const std::vector<double> & reference) const;

        /*
        /// @brief get element (i, j) of a matrix stored in a vector
        /// @param i row index
        /// @param j column index
        /// @return element (i, j) of the computed solution
        /// @note non-const version of the operator()
        double operator()(size_t i, size_t j) {
            return sol[i * n + j];
        };

        /// @brief get element (i, j) of a matrix stored in a vector
        /// @param i row index
        /// @param j column index
        /// @return element (i, j) of the exact solution
        /// @note const version of the operator()
        const double operator()(size_t i, size_t j) const {
            return sol[i * n + j];
        };
        */

        /// @brief print the computed solution
        void printSolution() const;

    private:
        /// @brief exact solution of the equation
        std::vector<double> exact_sol;

        /// @brief computed solution
        std::vector<double> sol;
        
        /// @brief right-hand side of the equation
        std::vector<double> rhs;

        /// @brief boundary conditions
        /// @details topbc and bottombc should have size n
        /// @details leftbc and rightbc should have size n-2
        std::vector<double> topbc;
        std::vector<double> rightbc;
        std::vector<double> bottombc;
        std::vector<double> leftbc;
        
        /// @brief number of grid points
        size_t n;
        
        /// @brief number of iterations
        unsigned n_iter = 0;
        
        /// @brief maximum number of iterations
        unsigned max_iter = 1000;
        
        /// @brief tolerance for convergence
        double tol = 1e-6;
};

#endif // SERIAL_SOLVER_HPP