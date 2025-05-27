set terminal png enhanced font 'Arial,12' size 800,600
set output 'test/plots/scalability.png'
set title 'Scalability Test'
set xlabel 'Number of Processes'
set ylabel 'Time (s)'
set grid
set key outside right
set logscale x 2
plot 'test/plots/scalability.dat' using 1:2 with linespoints title 'n=56', 'test/plots/scalability.dat' using 1:3 with linespoints title 'n=64'
