for Nb in 31 63
do
for Na in 0 1 2 3 4 5  
do
dir=SIM_NB_"$Nb"_NA_"$Na"
echo $dir
mkdir $dir
cd $dir 

cp ../brush.tcl .
sed -i "s/N\_B 31/N\_B $Nb/" brush.tcl
sed -i "s/N\_A 1/N\_A $Na/" brush.tcl
Espresso brush.tcl > output 2>&1 & wait %1;
process.sh > data.dat 2>&1

cd ..
done
done

#cd ..
#cd TEST_cshift
#cd CHARGED
#./simulate.sh
#cd ..
#cd UNCHARGED
#./simulate.sh
#cd ..

