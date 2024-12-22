set xlabel "Омега"
set ylabel "Итерации"

set grid

set xrange [0:2]
set yrange [0:24000]

plot "omega.csv" using 1:2 with linespoints

set terminal png
set output "omega.png"

replot

set terminal x11

pause -1

