#!/bin/bash

make

# Scalability test script using mpirun with different processor counts

echo ""
echo ""
echo "=================================================================="
echo "==== Running scalability test with different processor counts ===="
echo "=================================================================="

# Test with 1 processor
echo ""
echo "========================"
echo "Testing with 1 processor"
echo "========================"
mpirun -np 1 ./main

# Test with 2 processors
echo ""
echo "========================="
echo "Testing with 2 processors"
echo "========================="
mpirun -np 2 ./main

# Test with 4 processors
echo ""
echo "========================="
echo "Testing with 4 processors"
echo "========================="
mpirun -np 4 ./main

# Test with 8 processors
echo ""
echo "========================="
echo "Testing with 8 processors"
echo "========================="
mpirun -np 8 ./main

echo ""
echo ""
echo "Scalability test completed."

gnuplot test/plots/l2error_vs_h.gp test/plots/l2error_vs_n.gp test/plots/timing_vs_h.gp test/plots/timing_vs_n.gp test/plots/scalability.gp

echo "Plots generated successfully."