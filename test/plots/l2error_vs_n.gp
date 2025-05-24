set terminal png enhanced font 'Arial,12' size 800,600
set output 'test/plots/l2error_vs_n.png'
set title 'L2 Error vs Grid Size (n)'
set xlabel 'n'
set ylabel 'L2 Error'
set grid
set key outside right
set logscale x 2
set logscale y 2
plot 'test/plots/l2error_vs_n.dat' using 1:2 with linespoints title 'L2 Error'
