Nb=31
Na=0
for eps in 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0
do
dir=SIM_NB_"$Nb"_NA_"$Na"_eps_"$eps"
cd $dir 

#process.sh > data.dat 2>&1
x=$(process.sh | head -n 2 | tail -n 1 | awk '{print $3}')
dx=$(process.sh | head -n 2 | tail -n 1 | awk '{print $5}')
echo "$eps $x $dx"
cd ..
done
