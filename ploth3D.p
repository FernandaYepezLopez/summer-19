set palette defined (0 'purple', 0.25 'blue', .5 'green', .75 'yellow', 1.0 'red')
set logscale cb
set cbrange[.0001:1]
set title 'GPD DGLAP H'
set xlabel 'X'
set ylabel 'Xi'
set zlabel 't'
set zrange[0:-2]
splot "outputh3D.dat" u 1:2:3:4 with points pointtype 7 ps 1 palette
