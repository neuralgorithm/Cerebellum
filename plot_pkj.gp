set size square

set xtics 0.5
set ytics 20

g1(x) = d1*cos(2*3.14*(x+e1)/2.0)+f1
g2(x) = d2*cos(2*3.14*(x+e2)/2.0)+f2
g3(x) = d3*cos(2*3.14*(x+e3)/2.0)+f3
g4(x) = d4*cos(2*3.14*(x+e4)/2.0)+f4

fit g1(x) '10_10.dat' via d1,e1,f1
fit g2(x) '100_10.dat' via d2,e2,f2
fit g3(x) '200_10.dat' via d3,e3,f3
fit g4(x) '300_10.dat' via d4,e4,f4

set terminal postscript eps
set output 'pkj_10.eps'
plot  [0:2][0:100] '10_10.dat' t '' w p, '100_10.dat' t '' w p, '200_10.dat' t '' w p, '300_10.dat' t '' w p, g1(x) t '1st' w l, g2(x) t '100th' w l, g3(x) t '200th' w l, g4(x) t '300th' w l
set terminal x11
replot
