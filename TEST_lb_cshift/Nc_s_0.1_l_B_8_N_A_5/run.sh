#!/bin/bash
### DIRECTORIES
HOME=/storage/praha1/home/kubeckaj
UTILITIES=$HOME/Utilities
WORKDIR=/storage/praha1/home/kubeckaj/BRUSHES/TEST_lb_cshift/Nc_s_0.1_l_B_8_N_A_5
### FILES OR EXECUTABLES
ESPRESSO=$UTILITIES/espresso-3.3.1/build_P3M/Espresso
WORKFILE=$WORKDIR/brush.tcl
COMPILEFILE=$UTILITIES/compileESPRESSO.sh


### RUN
$($COMPILEFILE)
cd $WORKDIR
mkdir RESULTS STRUCTURES FUNCTIONS CORRELATIONS
#$ESPRESSO $WORKFILE > output 2>&1
rm fft*
cd RESULTS
python $UTILITIES/ShapeAnisotropyJK anisotropy_backbone.txt
python $UTILITIES/ShapeAnisotropyJK anisotropy_all.txt
for i in `ls`; do python average_t_X_Y_Z.py $i; done
cd ..
cd FUNCTIONS
for i in `ls`; do python average_ttt_XXX.py $i; done
cd ..
cd CORRELATIONS
for i in `ls`; do python average_t_sl_XXX.py $i; done
cd ..
