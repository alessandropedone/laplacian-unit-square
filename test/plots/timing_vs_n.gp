set terminal png enhanced font 'Arial,12' size 800,600
set output 'test/plots/timing_vs_n.png'
set title 'Timing vs Grid Size (n)'
set xlabel 'n'
set ylabel 'Time (s)'
set grid
set key outside right
set logscale x 2
set logscale y 2
plot 'test/plots/timing_vs_n.dat' using 1:2 with linespoints title 'Serial', 'test/plots/timing_vs_n.dat' using 1:3 with linespoints title 'OMP', 'test/plots/timing_vs_n.dat' using 1:4 with linespoints title 'MPI', 'test/plots/timing_vs_n.dat' using 1:5 with linespoints title 'Hybrid'
