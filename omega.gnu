set xlabel "Omega"
set ylabel "Iterations"

set grid

set xrange [1:2]
set yrange [0:8000]

plot "omega.csv" using 1:2 with linespoints

set terminal png
set output "omega.png"

replot

set terminal x11

pause -1

