[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/bOfolMCC)
# A matrix–free parallel solver for the Laplace equation

## Roadmap

### Extras
- use getpot to read data
- (Neumann and Robin boundary conditions)[https://chatgpt.com/share/68322bf2-c650-8006-abfa-ca0234cb86ef]

### ?
- Schwarz iteration type (change only the local solver)

### README
- results discussion
- instruction to reproduce the scalability and grid size tests, just run
    ```bash
    ./test.sh
    ```
- explain generated plots
- simpler python script to visdualize results
- hardware comparison instructions
- instruction for compilation (OPENMP = 1 option)
- required package (muparserx), run if on Debian-based Linux distribution
    ```bash
    sudo apt-get update
    sudo apt-get install libmuparserx-dev
    ```



