#!/bin/bash

source setup.sh ew
source slabname.sh

echo "%%%%%%%%%%% STEP A12 - show the density profiles of all simulations %%%%%%%%%%%"
if [ ! -e $SLAB.cp ] ; then
  echo "ERROR! File $SLAB.cp is missing."
  echo "Maybe the simulation is still running - run Step 11 to interrupt."
  sleep 3
  exit
fi

rm *~

nice showcp -a -b10 $SLAB.cp

cat <<EOF
This task is designed for a group of students working in parallel in
subdirectories, each student running the simulations at different temperature.  This scrip shows all obtained convergence profiles
(density vs. time)."
EOF

plot :0:4 slab_*K.cpa ../*/slab_*K.cpa $SLAB.cpa:0:4:=1

echo
echo "                       \\\\\\\\\\V/////"
echo "                      /  _    _  \\"
echo "                     |   O    O   |"
echo "                    @|     ..     |@"
echo "                      \\   \\__/   /"
echo "                       \\________/"
echo
