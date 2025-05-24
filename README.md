[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/bOfolMCC)
# A matrix–free parallel solver for the Laplace equation

## Roadmap

### Jacobi iteration method

Test:
- L2 norm of the error as function of grid size (grafichino)
- f(x) = 8π2 sin(2πx) sin(2πy), with exact sol u(x, y) = sin(2πx) sin(2πy)
- serial vs parallel as grid size increases (time using chrono)
- test folder: scalability test (bash script) (1,2,4 cores), data folder to collect results

### Extras
- pybind11 to plot [results](https://chatgpt.com/share/682df901-b6cc-8006-9be6-a300f33212a8)
- use getpot to read data
- (Neumann and Robin boundary conditions)[https://chatgpt.com/share/68322bf2-c650-8006-abfa-ca0234cb86ef]
- Schwarz iteration type (change only the local solver)

### README
- results discussion
- instruction to reproduce the scalability test
- hardware comparison instructions
- instruction for compilation (OPENMP = 1 option)