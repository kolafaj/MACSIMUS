#!/bin/bash
if [ $# -lt 1 ]; then 
  echo "Viscosity analysis. Call by:"
  echo "  visc.sh SIMNAME [LINE]"
  echo "LINE=0 = last line (of .run for final results)"
  echo "NB: the error is unreliable, on average should be multiplied by 1.3"
  exit
fi

staprt -t2 \
  -n'Ptyz [Pa]' \
  -n'Ptzx [Pa]' \
  -n'Ptxy [Pa]' \
  -n'Ptxx-tr [Pa]' \
  -n'Ptyy-tr [Pa]' \
  -n'Ptzz-tr [Pa]' \
  $1 

runint -m2 -o < $1.PtyzPa.cov > $1.PtyzPa.run
runint -m2 -o < $1.PtzxPa.cov > $1.PtzxPa.run
runint -m2 -o < $1.PtxyPa.cov > $1.PtxyPa.run
runint -m2 -o < $1.Ptxx-trPa.cov | tabproc A 'B*.75' > $1.Ptxx-trPa.run
runint -m2 -o < $1.Ptyy-trPa.cov | tabproc A 'B*.75' > $1.Ptyy-trPa.run
runint -m2 -o < $1.Ptzz-trPa.cov | tabproc A 'B*.75' > $1.Ptzz-trPa.run

mergetab \
  $1.PtyzPa.run:1:2 \
  $1.PtzxPa.run:1=1:2 \
  $1.PtxyPa.run:1=1:2 \
  $1.Ptxx-trPa.run:1=1:2 \
  $1.Ptyy-trPa.run:1=1:2 \
  $1.Ptzz-trPa.run:1=1:2 | \
  tabproc A 'B/5+C/5+D/5+E/7.5+F/7.5+G/7.5' \
    '\(B^2/5+C^2/5+D^2/5+E^2/7.5+F^2/7.5+G^2/7.5-(B/5+C/5+D/5+E/7.5+F/7.5+G/7.5)^2)/2' > $1.Pav.run
# /2 is for /(5-1), where 5 is the # of indep. measurements
# NB: the error is underestimated, should be multiplied by 1.3

plot $1.PtyzPa.run:1:2:-2 \
     $1.PtzxPa.run:1:2:-2 \
     $1.PtxyPa.run:1:2:-2 \
     $1.Ptxx-trPa.run:1:2:-3 \
     $1.Ptyy-trPa.run:1:2:-3 \
     $1.Ptzz-trPa.run:1:2:-3 \
     $1.Pav.run:1:2:-1:3

if [ "$2" == "" ] ; then L=0; else L=$2 ; fi

if [ "$L" == "0" ] ; then
  tail -q -n1 $1.PtyzPa.run $1.PtzxPa.run $1.PtxyPa.run | sumetc -a | field 2 > aux
  tail -q -n1 $1.Ptxx-trPa.run $1.Ptyy-trPa.run $1.Ptzz-trPa.run | sumetc -a | field 2 >> aux
  tail -q -n1 $1.Pav.run | tabproc B C > auxx
else
  cat $1.PtyzPa.run $1.PtzxPa.run $1.PtxyPa.run | line $L | sumetc -a | field 2 > aux
  cat $1.Ptxx-trPa.run $1.Ptyy-trPa.run $1.Ptzz-trPa.run | line $L | sumetc -a | field 2 >> aux
  line $L < $1.Pav.run | tabproc B C > auxx
fi
oneline < aux | tabproc "0.6*A+0.4*B" >> aux
cat auxx >> aux
echo "# line=$L off-diag traceless FINAL FINAL stderr" >> aux
oneline aux
