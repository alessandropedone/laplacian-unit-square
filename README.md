[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/bOfolMCC)
# A matrix–free parallel solver for the Laplace equation

## Refinements
- manage test for n = 2^7, 2^8 (requested)
    - per n=8 raggiunge max_iter in alcuni processori ma ottiene lo stesso il risultato finale
- review muparserx interface for speed
- review documentation and comments
- non convergence warnings
- negative times
- review gaphs


### README
- strucutre of the repo and doxygen documentation
- clone with submodules
- results discussion
- instruction to reproduce the scalability and grid size tests, just run
    ```bash
    ./test.sh
    ```
- explain generated plots
- simpler python script to visdualize results
- hardware comparison instructions
- instruction for compilation (OPENMP = 1 option)
- required package (muparserx and ltbb), run if on Debian-based Linux distribution
    ```bash
    sudo apt-get update
    sudo apt-get install libmuparserx-dev
    sudo apt-get install libtbb-dev
    ```
- flag --use-datafile to use dataset file information
- if you don't have 8 cores on your machine modify test.sh



