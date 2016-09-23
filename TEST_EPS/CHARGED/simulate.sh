Nb=31
Na=0
for eps in 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0
do
dir=SIM_NB_"$Nb"_NA_"$Na"_eps_"$eps"
echo $dir
mkdir $dir
cd $dir 

cp ../brush.tcl .
cp ../process.sh .
sed -i "s/N\_B 31/N\_B $Nb/" brush.tcl
sed -i "s/N\_A 1/N\_A $Na/" brush.tcl
sed -i "s/ eps 1.0/ eps $eps/" brush.tcl
Espresso brush.tcl > output 2>&1 & wait %1;
chmod +x process.sh
./process.sh > data.dat 2>&1

cd ..
done
