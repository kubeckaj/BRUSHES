reset 
set terminal png enhanced truecolor
set output "out.png"

FROM=500
f(x) = mean_y
fit [FROM:] f(x) 'results.txt' u 1:3 every 100 via mean_y
stddev_y = sqrt(FIT_WSSR / (FIT_NDF + 1))

f2(x) = (x > FROM ) ? mean_y : 1/0 
f2a(x) = (x > FROM ) ? mean_y+stddev_y : 1/0
f2b(x) = (x > FROM ) ? mean_y-stddev_y : 1/0

set arrow from FROM,0 to FROM,30 nohead lw 2.0 front
set yrange [0:30]
#set xrange [0:100]


p 'results.txt' u 1:3 t 'end-to-end distance'\
, f2(x) lw 3.0 t 'mean'\
, f2a(x) lt 0 lw 2.0 lc -1 t 'std.dev.'\
, f2b(x) lt 0 lw 2.0 lc -1 t ''



set out
exit
print mean_y
print stddev_y
