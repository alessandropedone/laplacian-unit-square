set terminal png enhanced font 'Arial,12' size 800,600
set output 'test/plots/l2error_vs_h.png'
set title 'L2 Error vs Grid Spacing (h)'
set xlabel 'h = 1/(n-1)'
set ylabel 'L2 Error'
set grid
set key outside right
set logscale x 2
set logscale y 2
plot 'test/plots/l2error_vs_h.dat' using 1:2 with linespoints title 'L2 Error'
