for Na in 0 1 3 5
do
Nb=31
dir=SIM_NB_"$Nb"_NA_"$Na"
cd $dir 
x=$(process.sh | head -n $1 | tail -n 1 | awk '{print $3}')
dx=$(process.sh | head -n $1 | tail -n 1 | awk '{print $5}')
echo "$Na $x $dx" >> ../dat1.dat
cd ..

Nb=63
dir=SIM_NB_"$Nb"_NA_"$Na"
cd $dir
x=$(process.sh | head -n $1 | tail -n 1 | awk '{print $3}')
dx=$(process.sh | head -n $1 | tail -n 1 | awk '{print $5}')
echo "$Na $x $dx" >> ../dat2.dat
cd ..
done
gnuplot G_yIxIdx.gp
gnome-open out.png
rm dat1.dat dat2.dat
