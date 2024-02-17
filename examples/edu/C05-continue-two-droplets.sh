#!/bin/bash

source setup.sh ce
source Cdata.sh

echo "%%%%%%%%%%%%%% STEP C05 - continue simulation of two droplets %%%%%%%%%%%%%%"

echo
echo "Use this to continue the two-droplet simulation"
echo "if it has not reached the final equilibrated state."
echo

if [ ! -e drop2.cfg ] ; then
  echo "ERROR! File drop2.cfg is missing. Re-run from STEP 03."
  echo "Type [Enter]"
  read
  exit
fi

cat > drop2.get <<EOF
thermostat=0
init="append"
no=n*12;
EOF

echo ">>>>>> Simulation parameters:"
cat drop2.get

echo
echo "The simulation of droplet coalescence will be performed in the background."
echo "It is possible to watch the running simulation (C04-watch-coalescence.sh)"
echo "[Enter]=continue, [Ctrl-C]=interrupt"
read

NAME=`basename $PWD`
$COOK spce drop2 &

echo
echo "Wait a minute or two and start Step C04"

echo "[Enter]=continue"
read
