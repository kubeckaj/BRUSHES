reset 
set terminal png enhanced truecolor
set output "out.png"

set xrange [-0.1:0.35]

p 'o' 	u 1:2:3 w yerr	lc 1 t 'file1'\
, 'o' 	u 1:2 w l	lc 1 t 'file1'\
, 'oo' 	u 1:2:3 w yerr	lc 2 t ''\
, 'oo' 	u 1:2 w l	lc 2 t ''\

set out
exit
