#!/bin/bash

source setup.sh ew

echo '%%%%%%%%%%% STEP A07 - simulation of the melt %%%%%%%%%%'

if [ ! -e melt.cfg ] ; then
  echo
  echo "ERROR! File melt.cfg is missing, repeat Step A06"
  sleep 3
  exit
fi

if [ -e melt.loc ] ; then
  echo "ERROR! The simulation of the same name is already running. I will try to interrupt it."
  int.sh melt
  echo "Finished! Type [Enter] to continue."
  read
  exit
fi

echo
echo "Now the molten salt simulation is about to be started."
echo "[Enter]=continue, [Ctrl-C]=interrupt"
read

cat > melt.get <<EOF
init="start" dt.plb=1 no=300;
EOF

$COOK nacl melt

echo "[Enter]=continue"
read
