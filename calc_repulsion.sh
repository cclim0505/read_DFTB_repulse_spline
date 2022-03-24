#!/bin/sh

echo "Reading $1 to calculate DFTB spline repulsion."

out_atomic="out_repulse_atomic.dat"
out_ev="out_repulse_eV.dat"

./a.out $1 > $out_atomic

awk '{printf"%4.7f \t\t %4.7f \n", $1*0.529177249, $2*27.211396641}' $out_atomic > $out_ev

