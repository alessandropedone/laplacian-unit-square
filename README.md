[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/bOfolMCC)
# A matrix–free parallel solver for the Laplace equation

## TBD
- review documentation and comments
- readme
    - non-convergin Scharz method
    - test all commands

## Setup and test
You can run the following command to clone the repository:
```bash
git clone --recurse-submodules git@github.com:PACS-24-25/challenge3-male.git
```
Then to compile and run all the tests you can execute:
```bash
./test.sh
```

## Required packages
It's required to link with `muparserx` library, which, for instance, if you're on Debian/Ubuntu/...  you can install on your machine with the following command:
```bash
sudo apt-get update
sudo apt-get install libmuparserx-dev
```

## Structure of the repository
```bash
challenge3-male
├── Challenge24-25-3.pdf
├── LICENSE
├── Makefile
├── README.md
├── data.txt
├── docs
│   ├── Doxyfile
│   ├── ale_hw.info
│   ├── compare_hw.sh
│   ├── generate_hw_info.sh
│   ├── html
│   └── marta_hw.info
├── include
│   ├── GetPot
│   ├── core
│   ├── eigen
│   ├── muparser_interface.hpp
│   ├── plot.hpp
│   ├── simulation_parameters.hpp
│   └── vtk.hpp
├── src
│   ├── main.cpp
│   └── solver.cpp
├── test
│   ├── data
│   ├── plot.py
│   └── plots
└── test.sh
```

### Docs
In `docs` folder you can find the documentation generated with doxygen in `html` format.

### Hardware information and comparison
Again in `docs` folder you can find the follwing files:
- `ale_hw.info` which contains the hardware information of [@alessandropedone](https://github.com/alessandropedone)
- `marta_hw.info` which contains the hardware information of [@martapignatelli](https://github.com/martapignatelli)
- `generate_hw_info.sh` which is an executable that you can run in the terminal to create a file with the information of your machine in the same format as the aforementioned files
- `compare_hw.sh` that allows to print a comparison between to hardware specifying the name of the files
For instance, you could run:
```bash
cd docs 
./generate_hw_info.sh
./compare_hw.sh ale_hw.info your_machine_name_hw.info
```

### Implementation
We implemented five methods to solve the problem, which are members of the class named Solver:
```cpp
/// @brief implement Jacobi iterative solver for the Laplace equation without parallelism
void solve_jacobi_serial();

/// @brief implement Jacobi iterative solver for the Laplace equation with OPENMP
void solve_jacobi_omp();

/// @brief implement Jacobi iterative solver for the Laplace equation with MPI
void solve_jacobi_mpi();

/// @brief implement Jacobi iterative solver for the Laplace equation with MPI and OpenMP
void solve_jacobi_hybrid();

/// @brief Schwarz implementation: locally, the equation is solved using Eigen LDLT decomposition
void solve_direct_mpi();
```
We chose to avoid using template programming since it would more complicated putting several `if constexpr` instead of simply structuring in a different way the code in separate member functions.

### Salability test
We performed a small scalability test with 1, 2 and 4 processors. \
The results can be obtained by running the command specified in the first section (timings are printed in the terminal).

### Grid size variation
We also made the grid size vary between 8 and 64, and we avoided going beyond this threshold because the execution took too long and results can be already observed with this choice of grid sizes.

## Flags
It's possible to disable the compilation with OPENMP by running
```bash
make OPENMP=0
```
There are two possibilities to run the code:
1. if you use the following command you just run the code on the example chosen by us,
    ```bash
    mpirun -np j ./main
    ```
2. if you run the same command but with `--use-datafile` flag the tests run with the data specified within `data.txt` (be careful: the code runs a lot slower because of the overhead of the interface).
    ```bash
    mpirun -np j ./main --use-datafile
    ```

## Results

In `test/data` folder you can find `.csv` files with saved timings from the last execution of the test and some saved solution in `.vtk` format. The results we obtained from running the test on our machine are already included in the repository. To view them, simply clone the repository without running the test again on yuor machine.

In `test/plot` you can find 5 plots in `.png` format, which are generated using gnuplot from `.csv` in `test/data` folder, of the following quantities:
- L2 error (one plot with n and one with h),
- timings of the case with two processors, stored in `results_2.csv` (one plot with n and one with h),
- scalability test for the hybrid method in the case of $n=56$ and $n=64$.

We also kept a python script, as an alternative. Only pandas and matplotlib libraries are required to run the python code.

The direct solver doesn't converge with 4 processors and finer mesh, and this is an intrinsic problem of Schwarz method, since... TBD!!! \
It's possible to observed the solution in the non convergence case of the direct solver in `test/data/solution_4_n_64.vtk`.


