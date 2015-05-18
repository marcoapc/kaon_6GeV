#!/bin/csh -f
set num = `echo $1 | cut -c7- | cut -d"_" -f1`
echo $num
grep -i "normfac " $1  | awk '{printf "%s\n",$3}' >! ./normfac_$num.dat
grep -i "charge " $1  | awk '{printf "%s\n",$3}' >! ./charge_$num.dat
grep -i "Ebeam " $1 | awk '{printf "%s\n",$3}' >! ./Ebeam_$num.dat
grep -i "angle " $1 | awk '{printf "%s\n",$3}' >! ./eTheta_$num.dat
grep -i "angle " $1 | awk '{printf "%s\n",$4}' >! ./pTheta_$num.dat
grep -i "momentum " $1 | awk '{printf "%s\n",$3}' >! ./eP_$num.dat
grep -i "momentum " $1 | awk '{printf "%s\n",$4}' >! ./pP_$num.dat
