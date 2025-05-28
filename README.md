[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/bOfolMCC)
# A matrix–free parallel solver for the Laplace equation

## Refinements
- review muparserx interface for speed
- review documentation and comments

### README
- structure of the repo and doxygen documentation
- clone with submodules
- results discussion
- instruction to reproduce the scalability and grid size tests, just run
    ```bash
    ./test.sh
    ```
- explain generated plots
- simpler python script to visualize results
- hardware comparison instructions
- instruction for compilation (OPENMP = 1 option)
- required package (muparserx and ltbb), run if on Debian-based Linux distribution
    ```bash
    sudo apt-get update
    sudo apt-get install libmuparserx-dev
    sudo apt-get install libtbb-dev
    ```
- flag --use-datafile to use dataset file information
- explain why direct solver doesn't converge with 4 processors and finer mesh (you can see the non converged solution in test/data/solution_4_n_64.vtk)
- too slow for grid size over 64

Plots:
- timing of result_2 (2 processors are plotted)
- scalability test for the hybrid method


