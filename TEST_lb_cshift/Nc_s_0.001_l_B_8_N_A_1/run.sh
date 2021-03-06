#!/bin/bash
### DIRECTORIES
HOME=/storage/praha1/home/kubeckaj
UTILITIES=$HOME/Utilities
WORKDIR=/storage/praha1/home/kubeckaj/BRUSHES/TEST_lb_cshift/Nc_s_0.001_l_B_8_N_A_1
### FILES OR EXECUTABLES
ESPRESSO=$UTILITIES/espresso-3.3.1/build_P3M/Espresso
WORKFILE=$WORKDIR/brush.tcl
COMPILEFILE=$UTILITIES/compileESPRESSO.sh


### RUN
$($COMPILEFILE)
cd $WORKDIR
mkdir RESULTS STRUCTURES FUNCTIONS CORRELATIONS
$ESPRESSO $WORKFILE > output 2>&1
rm fft*
