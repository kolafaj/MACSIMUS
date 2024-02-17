#!/bin/bash

source setup.sh ce
source Cdata.sh

echo "%%%%%%%%%%% STEP C04 - watching the coalescence %%%%%%%%%%%"

if [ ! -s drop2.plb ] ; then
  echo "ERROR! File slab.plb not (yet) available. Wait and re-try."
  echo "[Enter]=continue"
  read
  exit
fi

export MOLCFG=1
molcfg -$N:a SPCE -$N:b SPCE drop2

nice show drop2 -s150 -Ca -Yb '-I%i'

if [ -e drop2.loc ] ; then
  if ./yes.sh "The simulation is still running. Interrupt (Y/n)?" y ; then
    int.sh drop2
  fi
fi

    echo
    echo "                       \\\\\\\\\\V/////"
    echo "                      /  _    _  \\"
    echo "                     |   O    O   |"
    echo "                    @|     ..     |@"
    echo "                      \\   \\__/   /"
    echo "                       \\________/"
    echo
rm *~
sleep 1
