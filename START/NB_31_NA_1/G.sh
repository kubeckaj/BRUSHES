#!/bin/bash

echo $1
gnuplot $1
gnome-open out.png
rm fit.log
