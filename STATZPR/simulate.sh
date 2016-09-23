for Na in 0 1 3 5
do
dir=SIM_NA_$Na
echo $dir
mkdir $dir
cd $dir 

cp ../brush.tcl .
cp ../process.sh .
sed -i "s/N\_A 1/N\_A $Na/" brush.tcl
Espresso brush.tcl > output 2>&1 & wait %1;
pids="$!"
wait $pids
chmod +x process.sh
./process.sh data.dat

cd ..
done
