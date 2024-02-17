#!/bin/bash

source setup.sh ew

echo "%%%%%%%%%%% STEP A05 - watching the results of Step A04 %%%%%%%%%%%"

if [ ! -e box.def ] ; then
  echo "ERROR! The expected files of Step A04 are missing. Repeat Step 04."
  sleep 3
  exit
fi

rm *~

./res.sh cryst rho
