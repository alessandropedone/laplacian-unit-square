/// @file jacobi_serial_solver.hpp
/// @brief Header file for the JacobiSerialSolver class
/// @details This file contains the declaration of the JacobiSerialSolver class,
///          which implements an iterative solver for a given equation.
///          The class provides methods to set the boundary conditions,
///          initial guess, exact solution, and right-hand side of the equation.
///          It also includes a method to print the computed solution.
/// @details The class is designed to be used in a serial computing environment.
///          It is not optimized for parallel computing and does not use any
///          parallel computing libraries or techniques.
/// @details The class is intended to be used as a comparison for more complex
///          solvers that may incorporate parallelism or other advanced features.
#ifndef JACOBI_SERIAL_SOLVER_HPP
#define JACOBI_SERIAL_SOLVER_HPP
#include <iostream>
#include <vector>
#include <functional>

class JacobiSerialSolver
{
public:
    JacobiSerialSolver() = delete;

    /// @brief constructor
    /// @param uex exact solution of the equation
    /// @param initial_guess initial guess for the solution
    /// @param f right-hand side of the equation
    /// @param top_bc top boundary condition
    /// @param right_bc right boundary condition
    /// @param bottom_bc bottom boundary condition
    /// @param left_bc left boundary condition
    /// @param n grid size
    /// @param max_iter maximum number of iterations
    /// @param tol tolerance for convergence
    JacobiSerialSolver(
        const std::vector<double> &initial_guess,
        std::function<double(double, double)> f,
        std::function<double(double, double)> top_bc,
        std::function<double(double, double)> right_bc,
        std::function<double(double, double)> bottom_bc,
        std::function<double(double, double)> left_bc,
        size_t n,
        unsigned max_iter = 1000,
        double tol = 1e-6,
        std::function<double(double, double)> uex = nullptr)
        : iter(0),
          n(n),
          max_iter(max_iter),
          tol(tol),
          uex(uex),
          uh(initial_guess),
          f(f),
          top_bc(top_bc),
          right_bc(right_bc),
          bottom_bc(bottom_bc),
          left_bc(left_bc)
    {
    }

    /// @brief default destructor
    ~JacobiSerialSolver() = default;

    /// @brief implement Jacobi iterative solver for the Laplace equation
    void solve();

    // SETTERS

    /// @brief set grid size
    /// @param n grid size
    /// @details The grid size is used to define the number of grid points
    ///          and the size of the solution vector
    void set_n(size_t n)
    {
        this->n = n;
    };

    /// @brief set the number of max iterations
    /// @param max_iter maximum number of iterations
    /// @details The maximum number of iterations is used to limit the number
    ///          of iterations in the iterative solver
    void set_max_iter(size_t max_iter)
    {
        this->max_iter = max_iter;
    };

    /// @brief set the tolerance for convergence
    /// @param tol tolerance for convergence
    /// @details The tolerance is used to check the convergence of the iterative
    ///          solver. If the residual is less than the tolerance,
    ///          the solver is considered to have converged
    void set_tol(double tol)
    {
        this->tol = tol;
    };

    /// @brief set the exact solution of the equation
    /// @param uex exact solution of the equation
    /// @details The exact solution is used to compare the computed solution
    ///          and to compute the error of the computed solution
    void set_uex(const std::function<double(double, double)> &uex)
    {
        this->uex = uex;
    };

    /// @brief set the initial guess for the solution
    /// @param initial_guess initial guess for the solution
    /// @details The initial guess is used to initialize the solution vector
    ///          before the iterative solver starts
    void set_initial_guess(const std::vector<double> &initial_guess)
    {
        this->uh = initial_guess;
    }

    /// @brief set the exact solution of the equation
    /// @param exact_sol exact solution of the equation
    void set_exact_sol(std::function<double(double, double)> uex)
    {
        this->uex = uex;
    };

    /// @brief set the right-hand side of the equation
    /// @param rhs right-hand side of the equation
    void set_f(std::function<double(double, double)> f)
    {
        this->f = f;
    };

    /// @brief set the boundary conditions of the system
    /// @param top_bc top boundary condition
    /// @param right_bc right boundary condition
    /// @param bottom_bc bottom boundary condition
    /// @param left_bc left boundary condition
    /// @details The boundary conditions are defined as functions of two variables
    ///          and are used to set the values of the solution at the boundaries
    void set_bc(const std::function<double(double, double)> &top_bc,
                const std::function<double(double, double)> &right_bc,
                const std::function<double(double, double)> &bottom_bc,
                const std::function<double(double, double)> &left_bc)
    {
        this->top_bc = top_bc;
        this->right_bc = right_bc;
        this->bottom_bc = bottom_bc;
        this->left_bc = left_bc;
    };

    // GETTERS

    /// @brief get the number of iterations tracked during the solver
    /// @return number of iterations
    unsigned get_iter() const
    {
        return iter;
    }

    /// @brief get the computed solution
    /// @return computed solution
    const std::vector<double> get_uh() const
    {
        return uh;
    };

    /// @brief get the exact solution in vector form
    /// @return exact solution
    const std::vector<double> get_uex() const
    {
        std::vector<double> temp(n * n);
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                temp[i * n + j] = uex(i, j);
            }
        }
        return temp;
    };

private:
    /// @brief number of iterations tracked during the solver
    unsigned iter = 0;

    /// @brief number of grid points
    size_t n;

    /// @brief maximum number of iterations
    unsigned max_iter = 1000;

    /// @brief tolerance for convergence
    double tol = 1e-6;

    /// @brief Exact solution of the equation
    /// @details uex should be a function of two variables
    std::function<double(double, double)> uex;

    /// @brief computed approximate solution of the equation
    /// @details uh should have size n*n
    std::vector<double> uh;

    /// @brief force term of the equation
    std::function<double(double, double)> f;

    /// @brief top boundary condition
    std::function<double(double, double)> top_bc;

    /// @brief right boundary condition
    std::function<double(double, double)> right_bc;

    /// @brief bottom boundary condition
    std::function<double(double, double)> bottom_bc;

    /// @brief left boundary condition
    std::function<double(double, double)> left_bc;

    /// @brief compute the L2 norm of the errror between two solutions in vector form
    /// @param sol1 first solution vector
    /// @param sol2 second solution vector
    /// @return error between the two solutions
    double compute_error(const std::vector<double> &sol1, const std::vector<double> &sol2) const;

    /// @brief compute the L2 norm of the error between a solution in vector form and one in function form
    /// @param sol1 computed solution vector
    /// @param sol2 solution as std::function (for example uex)
    /// @return error between the two solutions
    double compute_error(const std::vector<double> &sol1, const std::function<double(double, double)> &sol2) const;

    /// @brief get element (i, j) of the computed solution
    /// @param i row index
    /// @param j column index
    /// @return element (i, j) of the computed solution
    double uh_at(size_t i, size_t j)
    {
        if (i < 0 || i >= n || j < 0 || j >= n)
        {
            throw std::out_of_range("Index out of range");
        }
        return uh[i * n + j];
    };

    /// @brief get element (i, j) of the exact solution, force term or boundary condition
    /// @param i row index
    /// @param j column index
    double fun_at(std::function<double(double, double)> fun, size_t i, size_t j) const
    {
        if (i < 0 || i >= n || j < 0 || j >= n)
        {
            throw std::out_of_range("Index out of range");
        }
        return fun(static_cast<double>(i) / (n - 1), static_cast<double>(j) / (n - 1));
    };
};

#endif // JACOBI_SERIAL_SOLVER_HPP