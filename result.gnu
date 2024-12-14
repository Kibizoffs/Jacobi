set title 'Solution'

set xlabel 'X'
set ylabel 'Y'
set zlabel 'Z'

set xrange [-0.1:1.1]
set yrange [-0.1:1.1]

set cblabel 'Z'
set palette defined \
( \
    0   'purple', \
    0.2 'blue', \
    0.4 'cyan', \
    0.6 'green', \
    0.8 'yellow', \
    1   'red' \
)
set colorbox

set output 'solution.png'
splot 'result.csv' every 8:8:8 using 1:2:3 with points pointtype 7 palette
unset output  # Clear output so gnuplot can move to the next plot

pause -1

set title 'Difference'

set xlabel 'X'
set ylabel 'Y'
set zlabel 'Z'

set xrange [-0.1:1.1]
set yrange [-0.1:1.1]

set cblabel 'Z'

set palette defined \
( \
    0.0 'yellow', \
    0.2 'orange', \
    1.0 'red' \
)
set colorbox

set output 'difference.png'
splot 'result.csv' every 8:8:8 using 1:2:4 with points pointtype 3 palette
unset output  # Clear output to stop writing to file

pause -1
