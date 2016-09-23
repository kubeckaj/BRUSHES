reset 
set terminal png enhanced truecolor
set output "out.png"

set xrange [-0.5:5.5]

p 'dat1.dat' 	u 1:2:3 w yerr lc 1 t '31'\
, 'dat2.dat' 	u 1:2:3 w yerr lc 2 t '63'\
, 'dat1.dat'    u 1:2   w l    lc 1 t ''\
, 'dat2.dat'    u 1:2   w l    lc 2 t ''

set out
exit
