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

splot 'result.csv' every 8:8:8 using 1:2:3 with points pointtype 7 palette

pause -1
