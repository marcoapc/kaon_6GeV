#!/bin/csh -f
set num = `echo $1 | cut -c4- | cut -d"." -f1`
echo $num
grep -i "bo2cur " gen$num.txt  | awk '{printf "%s\n",$3}' >! ./Scalars/charge_$num.dat
grep -i "bo2cur " gen$num.txt  | awk '{printf "%s\n",$2}' >! ./Scalars/current_$num.dat
grep -i "bocur " gen$num.txt  | awk '{printf "%s\n",$3}' >! ./Scalars/charge2_$num.dat
grep -i "bocur " gen$num.txt  | awk '{printf "%s\n",$2}' >! ./Scalars/current2_$num.dat
#grep -i "*HMS fid"  gen$num.txt  | awk '{printf "%s\n",$5}' >! ./Scalars/htrack_$num.dat
grep -i "E SING FID TRACK EFFIC" hms$num.txt  | awk '{printf "%s\n",$7}' >! ./Scalars/htrack_$num.dat
grep "*SOS fid" gen$num.txt  | awk '{printf "%s\n",$5}' >! ./Scalars/strack_$num.dat
grep HMSTRG gen$num.txt  | awk '{printf "%s\n",$3}' >! ./Scalars/trig_$num.dat
grep HMSPRE gen$num.txt  | awk '{printf "%s\n",$3}' >! ./Scalars/ptrig_$num.dat
grep  "hpre" gen$num.txt  | awk '{printf "%s\n",$5}' >! ./Scalars/ptrig2_$num.dat
grep SOSPRE gen$num.txt  | awk '{printf "%s\n",$3}' >! ./Scalars/sptrig_$num.dat
grep  "spre" gen$num.txt  | awk '{printf "%s\n",$5}' >! ./Scalars/sptrig2_$num.dat
grep hELCLEAN gen$num.txt  | awk '{printf "%s\n",$5}' >! ./Scalars/helclean_$num.dat
grep "hCER " gen$num.txt  | awk '{printf "%s\n",$5}' >! ./Scalars/hcer_$num.dat
grep "hCERHI" gen$num.txt  | awk '{printf "%s\n",$5}' >! ./Scalars/hcerhi_$num.dat
grep "hPRE50" gen$num.txt | awk '{printf "%s\n",$3}' >! ./Scalars/hpre50_$num.dat
grep "hPRE100" gen$num.txt | awk '{printf "%s\n",$3}' >! ./Scalars/hpre100_$num.dat
grep "hPRE150" gen$num.txt | awk '{printf "%s\n",$3}' >! ./Scalars/hpre150_$num.dat
grep "*PS1(programmed)" gen$num.txt | awk '{printf "%s\n",$2}' >! ./Scalars/ps1_$num.dat
grep "*PS2(programmed)" gen$num.txt | awk '{printf "%s\n",$2}' >! ./Scalars/ps2_$num.dat
grep "PS3: " gen$num.txt | awk '{printf "%s\n",$8}' >! ./Scalars/ps3_$num.dat
grep -i "Beam on time 2" gen$num.txt  | awk '{printf "%s\n",$6}' >! ./Scalars/botime_$num.dat
grep -i "Beam on time 1" gen$num.txt  | awk '{printf "%s\n",$6}' >! ./Scalars/botime2_$num.dat
grep "coin adc" gen$num.txt | awk '{printf "%s\n",$4}' >! ./Scalars/ctrig_$num.dat
#grep "cTRG2" gen$num.txt | awk '{printf "%s\n",$3}' >! ./Scalars/ctrig_$num.dat
grep "all adcgates" gen$num.txt | awk '{printf "%s\n",$3}' >! ./Scalars/alladc_$num.dat
grep "sos adcgates" gen$num.txt | awk '{printf "%s\n",$3}' >! ./Scalars/sosadc_$num.dat
grep -i "Beam on time 2" gen$num.txt  | awk '{printf "%s\n",$9}' >! ./Scalars/botime_$num.dat
# Marco
grep -i "hcer e- cut eff=" gen$num.txt | awk '{printf "%s\n",$5}' >! ./Scalars/hcer_eEff_$num.dat
grep -i "hcal e- cut eff=" gen$num.txt | awk '{printf "%s\n",$5}' >! ./Scalars/hcal_eEff_$num.dat
# grep -i "PS1(calc)" gen$num.txt | awk '{printf "%s\n",$3}' >! ./Scalars/PS1_$num.dat
grep -i "E_beam " gen$num.txt | awk '{printf "%s\n",$3}' >! ./Scalars/Ebeam_$num.dat
grep -i "P HMS " gen$num.txt | awk '{printf "%s\n",$4}' >! ./Scalars/pHMS_$num.dat
grep -i "Theta HMS " gen$num.txt | awk '{printf "%s\n",$4}' >! ./Scalars/thetaHMS_$num.dat
grep -i "P SOS " gen$num.txt | awk '{printf "%s\n",$4}' >! ./Scalars/pSOS_$num.dat
grep -i "Theta SOS " gen$num.txt | awk '{printf "%s\n",$4}' >! ./Scalars/thetaSOS_$num.dat
grep -i "Beam on time 2" gen$num | awk '{print $6}' >! ./Scalers/bot_$num.dat
