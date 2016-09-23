#!/bin/bash
if [ -z $1 ];then 
echo "HELP: Please use this input"
echo "      ./P_yIxI.sh [number of observable]";
exit;
fi
Nb=31
Na=$2
for cshift in 0.0001 0.05 0.1 0.15 0.2 0.25
do
dir=SIM_NB_"$Nb"_NA_"$Na"_cshift_"$cshift"
cd $dir 

#process.sh > data.dat 2>&1
x=$(process.sh | head -n $1 | tail -n 1 | awk '{print $3}')
dx=$(process.sh | head -n $1 | tail -n 1 | awk '{print $5}')
echo "$cshift $x $dx"
cd ..
done
