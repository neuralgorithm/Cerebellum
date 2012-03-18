set size square

set xtics 0.5
set ytics 20

f1(x) = a1*cos(2*3.14*(x+b1)/2.0)+c1
f2(x) = a2*cos(2*3.14*(x+b2)/2.0)+c2
f3(x) = a3*cos(2*3.14*(x+b3)/2.0)+c3
f4(x) = a4*cos(2*3.14*(x+b4)/2.0)+c4

fit f1(x) '10_16.dat' via a1,b1,c1
fit f2(x) '100_16.dat' via a2,b2,c2
fit f3(x) '200_16.dat' via a3,b3,c3
fit f4(x) '300_16.dat' via a4,b4,c4

plot [0:2][0:130] '10_16.dat' t '' w p, '100_16.dat' t '' w p, '200_16.dat' t '' w p, '300_16.dat' t '' w p, f1(x) t '1st' w l, f2(x) t '100th' w l, f3(x) t '200th' w l, f4(x) t '300th' w l
