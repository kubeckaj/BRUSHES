reset 
set terminal png enhanced truecolor
set output "out.png"

set xrange [-0.5:2.5]
p 'CHARGED/D_epsIreIrex.dat' 	u 1:2:3 w yerr t 'charged'\
, 'UNCHARGED/D_epsIreIrex.dat' 	u 1:2:3 w yerr t 'uncharged'
set out
exit
