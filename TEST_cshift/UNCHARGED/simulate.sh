Nb=31
for Na in 0 5
do
for cshift in 0.0001 0.05 0.1 0.15 0.2 0.25
do
dir=SIM_NB_"$Nb"_NA_"$Na"_cshift_"$cshift"
echo $dir
mkdir $dir
cd $dir 

cp ../brush.tcl .
sed -i "s/N\_B 31/N\_B $Nb/" brush.tcl
sed -i "s/N\_A 0/N\_A $Na/" brush.tcl
sed -i "s/Nc\_s 0.25/Nc\_s $cshift/" brush.tcl
Espresso brush.tcl > output 2>&1 & wait %1;
process.sh > data.dat 2>&1

cd ..
done
done
