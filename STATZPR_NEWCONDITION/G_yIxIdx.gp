reset 
set terminal png enhanced truecolor
set output "out.png"

set xrange [-0.5:5.5]

p 'dat1.dat' 	u 1:2:3 w yerr t '31'\
, 'dat2.dat' 	u 1:2:3 w yerr t '63'

set out
exit
