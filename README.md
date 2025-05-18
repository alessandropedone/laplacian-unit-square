[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/bOfolMCC)
# A matrix–free parallel solver for the Laplace equation

## Roadmap

### Jacobi iteration method

Data structure: std::vector (is it simple to export in vtk fromat?)

Serial solver (local solver)

Parallel solver
- local solver (parallelize with OpenMP)
- global solver (MPI, split nodes by rows, communication between processes)

Error, convergence (converged if all ranks satisfy the criterion) and max iterations (slow convergence!!)

Export the solution in vtk format

Test:
- f(x) = 8π2 sin(2πx) sin(2πy), with exact sol u(x, y) = sin(2πx) sin(2πy)
- serial vs parallel as grid size increases (time)
- L2 norm of the error as function of grid size (grafichino)
- test folder: scalability test (bash script) (1,2,4 cores), data folder to collect results
- It could be interesting also to report your hardware in a file hw.info: cat /proc/cpuinfo.

### Extras
- non homogeneous Dirichlet boundary conditions
- Neumann and Robin boundary conditions
- Schwarz iteration type (change only the local solver)

### README
- results discussion
- instruction to reproduce the scalability test